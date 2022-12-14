import time

#import h5py
#import hdbscan
import numpy as np
import torch
import SimpleITK as sitk
import os
import scipy.ndimage
import skimage.transform

from sklearn.cluster import MeanShift

from pytorch3dunet.datasets.utils import SliceBuilder
from pytorch3dunet.unet3d.utils import get_logger
from pytorch3dunet.unet3d.utils import remove_halo

logger = get_logger('UNet3DPredictor')


class _AbstractPredictor:
    def __init__(self, model, loader, output_file, config, **kwargs):
        self.model = model
        self.loader = loader
        self.output_file = output_file
        self.config = config
        self.predictor_config = kwargs

    @staticmethod
    def _volume_shape(dataset):
        # TODO: support multiple internal datasets
        raw = dataset.raws[0]
        if raw.ndim == 3:
            return raw.shape
        else:
            return raw.shape[1:]

    @staticmethod
    def _get_output_dataset_names(number_of_datasets, prefix='predictions'):
        if number_of_datasets == 1:
            return [prefix]
        else:
            return [f'{prefix}{i}' for i in range(number_of_datasets)]

    def predict(self):
        raise NotImplementedError


class StandardPredictor(_AbstractPredictor):
    """
    Applies the model on the given dataset and saves the result in the `output_file` in the H5 format.
    Predictions from the network are kept in memory. If the results from the network don't fit in into RAM
    use `LazyPredictor` instead.

    The output dataset names inside the H5 is given by `des_dataset_name` config argument. If the argument is
    not present in the config 'predictions{n}' is used as a default dataset name, where `n` denotes the number
    of the output head from the network.

    Args:
        model (Unet3D): trained 3D UNet model used for prediction
        data_loader (torch.utils.data.DataLoader): input data loader
        output_file (str): path to the output H5 file
        config (dict): global config dict
    """

    def __init__(self, model, loader, output_file, config, **kwargs):
        super().__init__(model, loader, output_file, config, **kwargs)
        self.it = 0

    def predict(self):
        out_channels = self.config['model'].get('out_channels')
        if out_channels is None:
            out_channels = self.config['model']['dt_out_channels']

        prediction_channel = self.config.get('prediction_channel', None)
        if prediction_channel is not None:
            logger.info(f"Using only channel '{prediction_channel}' from the network output")

        device = self.config['device']
        output_heads = self.config['model'].get('output_heads', 1)

        logger.info(f'Running prediction on {len(self.loader)} batches...')

        # dimensionality of the the output predictions
        # volume_shape = self._volume_shape(self.loader.dataset)
        x, _ = self.loader.dataset[0]
        # volume_shape = (x.shape)[1:]
        volume_shape = (296, 296, 296)

        if prediction_channel is None:
            prediction_maps_shape = (out_channels,) + volume_shape
        else:
            # single channel prediction map
            prediction_maps_shape = (1,) + volume_shape

        logger.info(f'The shape of the output prediction maps (CDHW): {prediction_maps_shape}')

        patch_halo = self.predictor_config.get('patch_halo', (8, 8, 5))
        self._validate_halo(patch_halo, self.config['loaders']['test']['slice_builder'])
        logger.info(f'Using patch_halo: {patch_halo}')

        # create destination H5 file
        # h5_output_file = h5py.File(self.output_file, 'w')
        np_output_file = np.zeros(volume_shape)
        # allocate prediction and normalization arrays
        logger.info('Allocating prediction and normalization arrays...')
        prediction_maps, normalization_masks = self._allocate_prediction_maps(prediction_maps_shape,
                                                                              output_heads, np_output_file)

        # Sets the module in evaluation mode explicitly (necessary for batchnorm/dropout layers if present)
        self.model.eval()
        # Set the `testing=true` flag otherwise the final Softmax/Sigmoid won't be applied!
        self.model.testing = True
        # Run predictions on the entire input dataset
        with torch.no_grad():
            # for batch, indices in self.loader:
            for batch, indices in self.loader:
                # indices = (indices,)
                # send batch to device
                batch = batch.to(device)

                # forward pass
                # batch = torch.unsqueeze(batch, 0)
                predictions = self.model(batch)

                # x = predictions.clone().cpu()
                # np.save('/home/shkim/Libraries/pytorch-3dunet/datasets/JW/prediction' + str(indices[0]),
                #         x.numpy())
                # logger.info(f'Prediciton for {indices[0]} saved as {"datasets/JW/prediction" + str(indices[0])}')
                # wrap predictions into a list if there is only one output head from the network
                if output_heads == 1:
                    predictions = [predictions]

                # for each output head
                for prediction, prediction_map, normalization_mask in zip(predictions, prediction_maps,
                                                                          normalization_masks):

                    # convert to numpy array
                    prediction = prediction.cpu().numpy()

                    # for each batch sample
                    for pred, index in zip(prediction, indices):
                        # save patch index: (C,D,H,W)
                        if prediction_channel is None:
                            channel_slice = slice(0, out_channels)
                        else:
                            channel_slice = slice(0, 1)
                        # print(index)
                        index = (channel_slice,) + index

                        if prediction_channel is not None:
                            # use only the 'prediction_channel'
                            logger.info(f"Using channel '{prediction_channel}'...")
                            pred = np.expand_dims(pred[prediction_channel], axis=0)

                        logger.info(f'Saving predictions for slice:{index}...')

                        # remove halo in order to avoid block artifacts in the output probability maps
                        u_prediction, u_index = remove_halo(pred, index, volume_shape, patch_halo)
                        # accumulate probabilities into the output prediction array
                        prediction_map[u_index] += u_prediction
                        # count voxel visits for normalization
                        normalization_mask[u_index] += 1

        # save results to
                    self._save_results(prediction_maps, normalization_masks, output_heads, np_output_file, self.loader.dataset)
        # close the output H5 file
        # h5_output_file.close()

    def _allocate_prediction_maps(self, output_shape, output_heads, output_file):
        # initialize the output prediction arrays
        prediction_maps = [np.zeros(output_shape, dtype='float32') for _ in range(output_heads)]
        # initialize normalization mask in order to average out probabilities of overlapping patches
        normalization_masks = [np.zeros(output_shape, dtype='uint8') for _ in range(output_heads)]
        return prediction_maps, normalization_masks

    def _save_results(self, prediction_maps, normalization_masks, output_heads, output_file, dataset):
        def _slice_from_pad(pad):
            if pad == 0:
                return slice(None, None)
            else:
                return slice(pad, -pad)

        # save probability maps
        prediction_datasets = self._get_output_dataset_names(output_heads, prefix='predictions')
        for prediction_map, normalization_mask, prediction_dataset in zip(prediction_maps, normalization_masks,
                                                                          prediction_datasets):
            prediction_map = prediction_map / normalization_mask

            if dataset.mirror_padding is not None:
                z_s, y_s, x_s = [_slice_from_pad(p) for p in dataset.mirror_padding]

                logger.info(f'Dataset loaded with mirror padding: {dataset.mirror_padding}. Cropping before saving...')

                prediction_map = prediction_map[:, z_s, y_s, x_s]

            logger.info(f'Saving predictions to: {self.output_file}/{prediction_dataset}...')
            # output_file.create_dataset(prediction_dataset, data=prediction_map, compression="gzip")
        # np.save('/home/shkim/Libraries/pytorch-3dunet/prediction/'+prediction_dataset, prediction_map)
        patient = os.listdir(os.path.join(os.getcwd(), '../../io', 'test'))
        self._save_dicom(prediction_map, os.path.join(os.getcwd(), '../../io', 'mask', patient[self.it]))
        self.it += 1
        
    @staticmethod
    def _save_dicom(newArray, filepath):
        def _load_template():
            patient = filepath.split('/')[-1]
            dir = os.path.join(os.getcwd(), '../../io', 'test', patient)
            # dir = '/home/shkim/SpyderProjects/0611/'
            # dir = '/home/shkim/Libraries/pytorch-3dunet/datasets/JW/test/'
            # dir = '/home/shkim/Libraries/pytorch-3dunet/datasets/H/test/'
            # patient = os.listdir(dir)
            # dir = os.path.join(dir, patient[0])
            logger.info("The template DCM directory is " + dir)
            assert os.path.isdir(dir), 'Cannot find the template directory'
            reader = sitk.ImageSeriesReader()
            dicomFiles = reader.GetGDCMSeriesFileNames(dir)
            reader.SetFileNames(dicomFiles)
            reader.MetaDataDictionaryArrayUpdateOn()
            reader.LoadPrivateTagsOn()
            reader.Execute()
            imgShape=[len(reader.GetFileNames()),
                      int(reader.GetMetaData(0, "0028|0010")),
                      int(reader.GetMetaData(0, "0028|0011"))]
            # image = reader.Execute()
            # img3d = sitk.GetArrayFromImage(image)
            # img3d = img3d.transpose((1, 2, 0))
            logger.info(f"Original image shape is {imgShape}")
            return reader, imgShape

        if not os.path.isdir(filepath):
            os.mkdir(filepath)

        reader, imgShape = _load_template()
        logger.info("The template is loaded.")
        
        # spacing_ratio = np.array([296, 296, 296], dtype=np.float64) / imgShape
        spacing_ratio = np.array([1, 1, 1], dtype=np.float64)
        logger.info(f"Spacaing ratio is {spacing_ratio}")
        newArray = np.squeeze(newArray)
        newArray = skimage.transform.resize(newArray, imgShape, anti_aliasing=False)
        newArray[newArray > 0.5] = 1.0
        newArray[newArray <= 0.5] = 0.0
        paddedArray = newArray.astype(np.int16)
        # paddedArray = newArray
        newImage = sitk.GetImageFromArray(paddedArray)
        newImage = sitk.Cast(newImage, sitk.sitkInt8)

        # newImage.CopyInformation(oldImage)
        writer = sitk.ImageFileWriter()
        writer.KeepOriginalImageUIDOn()
        
        sp_x, sp_y = reader.GetMetaData(0, "0028|0030").split('\\')
        # sp_z = reader.GetMetaData(0, "0018|0050")
        _, _, z_0 = reader.GetMetaData(0, "0020|0032").split('\\')
        _, _, z_1 = reader.GetMetaData(1, "0020|0032").split('\\')
        sp_z = abs(float(z_0) - float(z_1))
    
        sp_z = float(sp_z) / spacing_ratio[0]
        sp_x = float(sp_x) / spacing_ratio[1]
        sp_y = float(sp_y) / spacing_ratio[2]

        modification_time = time.strftime("%H%M%S")
        modification_date = time.strftime("%Y%m%d")
        direction = newImage.GetDirection()
        series_tag_values = [(k, reader.GetMetaData(0, k)) for k in reader.GetMetaDataKeys(0)] + \
                             [("0008|0031", modification_time),
                             ("0008|0021", modification_date),
                             ("0028|0100", "8"),
                             ("0028|0101", "8"),
                             ("0028|0102", "7"),
                             ("0028|0103", "1"),
                             ("0028|0002", "1"),
                             ("0008|0008", "DERIVED\\SECONDARY"),
                             ("0020|000e", "1.2.826.0.1.3680043.2.1125." + modification_date + ".1" + modification_time),
                             ("0020|0037", '\\'.join(map(str, (direction[0], direction[3], direction[6], direction[1], direction[4], direction[7]))))]
    #                         ("0008|103e", reader.GetMetaData(0, "0008|103e") + " Processed-SimpleITK")]
    #    print(series_tag_values)
        
        logger.info(f'Saving mask into {filepath}')
        tags_to_skip = ['0010|0010', '0028|0030', '7fe0|0010', '7fe0|0000', '0028|1052',
                        '0028|1053', '0028|1054', '0010|4000', '0008|1030', '0010|1001',
                        '0008|0080', '0010|0040', '0008|1010']
        for i in range(newImage.GetDepth()):
            image_slice = newImage[:, :, i]
            # image_slice.CopyInformation(oldImage[:, :, i])
            for tag, value in series_tag_values:
                if (tag in tags_to_skip):
                    continue                
                if i == 0:
                    try:
                        logger.info(f'{tag} | {value}')
                    except:
                        continue
                image_slice.SetMetaData(tag, value)
            image_slice.SetMetaData("0008|0012", time.strftime("%Y%m%d"))
            image_slice.SetMetaData("0008|0013", time.strftime("%H%M%S"))
            image_slice.SetMetaData("0020|0032", reader.GetMetaData(i, "0020|0032"))
            # image_slice.SetMetaData('0020|0032', '\\'.join(map(str, [0, 0, i*sp_z])))
            image_slice.SetMetaData("0020|0013", str(i))
            image_slice.SetMetaData('0028|0030', '\\'.join(map(str, [sp_x, sp_y])))
            image_slice.SetSpacing([sp_x, sp_y])
            image_slice.SetMetaData("0018|0050", str(sp_z))
            writer.SetFileName(os.path.join(filepath, str(i).zfill(3)+'.dcm'))
            writer.Execute(image_slice)
        logger.info(f'Saved mask into {filepath}')
        
    @staticmethod
    def _save_dicom_legacy(newArray, filepath):
        def _load_template():
            patient = filepath.split('/')[-1]
            dir = os.path.join(os.getcwd(), '../../io', 'test', patient)
            # dir = '/home/shkim/SpyderProjects/0611/'
            # dir = '/home/shkim/Libraries/pytorch-3dunet/datasets/JW/test/'
            # dir = '/home/shkim/Libraries/pytorch-3dunet/datasets/H/test/'
            # patient = os.listdir(dir)
            # dir = os.path.join(dir, patient[0])
            logger.info("The template DCM directory is " + dir)
            assert os.path.isdir(dir), 'Cannot find the template directory'
            reader = sitk.ImageSeriesReader()
            dicomFiles = reader.GetGDCMSeriesFileNames(dir)
            reader.SetFileNames(dicomFiles)
            reader.MetaDataDictionaryArrayUpdateOn()
            reader.LoadPrivateTagsOn()
            image = reader.Execute()
            img3d = sitk.GetArrayFromImage(image)
            # img3d = img3d.transpose((1, 2, 0))
            return image, reader, img3d.shape

        if not os.path.isdir(filepath):
            os.mkdir(filepath)

        oldImage, reader, imgShape = _load_template()
        logger.info("The template is loaded.")
        newArray = np.squeeze(newArray)
        # newArray = scipy.ndimage.zoom(newArray, 2, order = 0)
        newArray = skimage.transform.resize(newArray, imgShape, anti_aliasing=False)
        # newArray = skimage.transform.rescale(newArray, 2, anti_aliasing = False)
        # newArray = skimage.transform.resize(newArray, (newArray.shape[0]*2, newArray.shape[1]*2, newArray.shape[2]*2), anti_aliasing=False)
        newArray[newArray > 0.5] = 1
        newArray[newArray <= 0.5] = 0
        # for _ in range(10):
        #     newArray = scipy.ndimage.binary_erosion(scipy.ndimage.binary_dilation(newArray))
        #     newArray = scipy.ndimage.binary_dilation(scipy.ndimage.binary_erosion(newArray))
        # newArray = scipy.ndimage.zoom(newArray, 2, order = 1)
        # paddedArray = np.pad(newArray.astype(np.int16), ((4, 4),)*3, 'constant', constant_values =  ((0, 0),)*3)
        paddedArray = newArray.astype(np.int16)

        # paddedArray = np.flip()
        # paddedArray = paddedArray.transpose((2, 0, 1))
        # paddedArray = paddedArray.transpose((1,2,0))
        # assert paddedArray.shape == (600,)*3, 'You idiot messed up with output dimension. Check unet3d/predictor.py --p.'
        newImage = sitk.GetImageFromArray(paddedArray)

        newImage.CopyInformation(oldImage)
        writer = sitk.ImageFileWriter()
        writer.KeepOriginalImageUIDOn()
        # tags_to_copy = [
        #                 "0010|0010",
        #                 "0010|0020",
        #                 "0010|0030",
        #                 "0020|000d",
        #                 "0020|0010",
        #                 "0008|0020",
        #                 "0008|0030",
        #                 "0008|0050",
        #                 "0008|0060",
        #                 ]
    #    tags_to_copy = list(reader.GetMetaDataKeys(0))

        modification_time = time.strftime("%H%M%S")
        modification_date = time.strftime("%Y%m%d")
        direction = newImage.GetDirection()
        series_tag_values = [(k, reader.GetMetaData(0, k)) for k in reader.GetMetaDataKeys(0)] + \
                            [("0008|0031", modification_time),
                             ("0008|0021", modification_date),
                             ("0008|0008", "DERIVED\\SECONDARY"),
                             ("0020|000e", "1.2.826.0.1.3680043.2.1125." + modification_date + ".1" + modification_time),
                             ("0020|0037", '\\'.join(map(str, (direction[0], direction[3], direction[6], direction[1], direction[4], direction[7]))))]
    #                         ("0008|103e", reader.GetMetaData(0, "0008|103e") + " Processed-SimpleITK")]
    #    print(series_tag_values)
        
        logger.info(f'Saving mask into {filepath}')
        for i in range(newImage.GetDepth()):
            image_slice = newImage[:, :, i]
            image_slice.CopyInformation(oldImage[:, :, i])
            for tag, value in series_tag_values:
                if i == 0:
                    logger.info(f'{tag}')
                if (tag == '0010|0010'):
                    continue
                image_slice.SetMetaData(tag, value)
            image_slice.SetMetaData("0008|0012", time.strftime("%Y%m%d"))
            image_slice.SetMetaData("0008|0013", time.strftime("%H%M%S"))
            image_slice.SetMetaData('0020|0032', '\\'.join(map(str, newImage.TransformIndexToPhysicalPoint((0, 0, i)))))
            image_slice.SetMetaData("0020|0013", str(i))
            writer.SetFileName(os.path.join(filepath, str(i).zfill(3) + '.dcm'))
            writer.Execute(image_slice)
        logger.info(f'Saved mask into {filepath}')

    @staticmethod
    def _validate_halo(patch_halo, slice_builder_config):
        patch = slice_builder_config['patch_shape']
        stride = slice_builder_config['stride_shape']

        patch_overlap = np.subtract(patch, stride)

        assert np.all(
            patch_overlap - patch_halo >= 0), f"Not enough patch overlap for stride: {stride} and halo: {patch_halo}"


class LazyPredictor(StandardPredictor):
    """
        Applies the model on the given dataset and saves the result in the `output_file` in the H5 format.
        Predicted patches are directly saved into the H5 and they won't be stored in memory. Since this predictor
        is slower than the `StandardPredictor` it should only be used when the predicted volume does not fit into RAM.

        The output dataset names inside the H5 is given by `des_dataset_name` config argument. If the argument is
        not present in the config 'predictions{n}' is used as a default dataset name, where `n` denotes the number
        of the output head from the network.

        Args:
            model (Unet3D): trained 3D UNet model used for prediction
            data_loader (torch.utils.data.DataLoader): input data loader
            output_file (str): path to the output H5 file
            config (dict): global config dict
        """

    def __init__(self, model, loader, output_file, config, **kwargs):
        super().__init__(model, loader, output_file, config, **kwargs)

    def _allocate_prediction_maps(self, output_shape, output_heads, output_file):
        # allocate datasets for probability maps
        prediction_datasets = self._get_output_dataset_names(output_heads, prefix='predictions')
        prediction_maps = [
            output_file.create_dataset(dataset_name, shape=output_shape, dtype='float32', chunks=True,
                                       compression='gzip')
            for dataset_name in prediction_datasets]

        # allocate datasets for normalization masks
        normalization_datasets = self._get_output_dataset_names(output_heads, prefix='normalization')
        normalization_masks = [
            output_file.create_dataset(dataset_name, shape=output_shape, dtype='uint8', chunks=True,
                                       compression='gzip')
            for dataset_name in normalization_datasets]

        return prediction_maps, normalization_masks

    def _save_results(self, prediction_maps, normalization_masks, output_heads, output_file, dataset):
        if dataset.mirror_padding:
            logger.warn(
                f'Mirror padding unsupported in LazyPredictor. Output predictions will be padded with pad_width: {dataset.pad_width}')

        prediction_datasets = self._get_output_dataset_names(output_heads, prefix='predictions')
        normalization_datasets = self._get_output_dataset_names(output_heads, prefix='normalization')

        # normalize the prediction_maps inside the H5
        for prediction_map, normalization_mask, prediction_dataset, normalization_dataset in zip(prediction_maps,
                                                                                                 normalization_masks,
                                                                                                 prediction_datasets,
                                                                                                 normalization_datasets):
            # split the volume into 4 parts and load each into the memory separately
            logger.info(f'Normalizing {prediction_dataset}...')

            z, y, x = prediction_map.shape[1:]
            # take slices which are 1/27 of the original volume
            patch_shape = (z // 3, y // 3, x // 3)
            for index in SliceBuilder._build_slices(prediction_map, patch_shape=patch_shape, stride_shape=patch_shape):
                logger.info(f'Normalizing slice: {index}')
                prediction_map[index] /= normalization_mask[index]
                # make sure to reset the slice that has been visited already in order to avoid 'double' normalization
                # when the patches overlap with each other
                normalization_mask[index] = 1

            logger.info(f'Deleting {normalization_dataset}...')
            del output_file[normalization_dataset]


class EmbeddingsPredictor(_AbstractPredictor):
    """
    Applies the embedding model on the given dataset and saves the result in the `output_file` in the H5 format.

    The resulting volume is the segmentation itself (not the embedding vectors) obtained by clustering embeddings
    with HDBSCAN or MeanShift algorithm patch by patch and then stitching the patches together.
    """

    def __init__(self, model, loader, output_file, config, clustering, iou_threshold=0.7, noise_label=-1, **kwargs):
        super().__init__(model, loader, output_file, config, **kwargs)

        self.iou_threshold = iou_threshold
        self.noise_label = noise_label
        self.clustering = clustering

        assert clustering in ['hdbscan', 'meanshift'], 'Only HDBSCAN and MeanShift are supported'
        logger.info(f'IoU threshold: {iou_threshold}')

        self.clustering_name = clustering
        self.clustering = self._get_clustering(clustering, kwargs)

    def predict(self):
        device = self.config['device']
        output_heads = self.config['model'].get('output_heads', 1)

        logger.info(f'Running prediction on {len(self.loader)} patches...')

        # dimensionality of the the output segmentation
        volume_shape = self._volume_shape(self.loader.dataset)

        logger.info(f'The shape of the output segmentation (DHW): {volume_shape}')

        logger.info('Allocating segmentation array...')
        # initialize the output prediction arrays
        output_segmentations = [np.zeros(volume_shape, dtype='int32') for _ in range(output_heads)]
        # initialize visited_voxels arrays
        visited_voxels_arrays = [np.zeros(volume_shape, dtype='uint8') for _ in range(output_heads)]

        # Sets the module in evaluation mode explicitly
        self.model.eval()
        self.model.testing = True
        # Run predictions on the entire input dataset
        with torch.no_grad():
            for batch, indices in self.loader:
                # logger.info(f'Predicting embeddings for slice:{index}')

                # send batch to device
                batch = batch.to(device)
                # forward pass
                embeddings = self.model(batch)

                # wrap predictions into a list if there is only one output head from the network
                if output_heads == 1:
                    embeddings = [embeddings]

                for prediction, output_segmentation, visited_voxels_array in zip(embeddings, output_segmentations,
                                                                                 visited_voxels_arrays):

                    # convert to numpy array
                    prediction = prediction.cpu().numpy()

                    # iterate sequentially because of the current simple stitching that we're using
                    for pred, index in zip(prediction, indices):
                        # convert embeddings to segmentation with hdbscan clustering
                        segmentation = self._embeddings_to_segmentation(pred)
                        # stitch patches
                        self._merge_segmentation(segmentation, index, output_segmentation, visited_voxels_array)

        # save results
        with h5py.File(self.output_file, 'w') as output_file:
            prediction_datasets = self._get_output_dataset_names(output_heads,
                                                                 prefix=f'segmentation/{self.clustering_name}')
            for output_segmentation, prediction_dataset in zip(output_segmentations, prediction_datasets):
                logger.info(f'Saving predictions to: {output_file}/{prediction_dataset}...')
                output_file.create_dataset(prediction_dataset, data=output_segmentation, compression="gzip")

    def _embeddings_to_segmentation(self, embeddings):
        """
        Cluster embeddings vectors with HDBSCAN and return the segmented volume.

        Args:
            embeddings (ndarray): 4D (CDHW) embeddings tensor
        Returns:
            3D (DHW) segmentation
        """
        # shape of the output segmentation
        output_shape = embeddings.shape[1:]
        # reshape (C, D, H, W) -> (C, D * H * W) and transpose -> (D * H * W, C)
        flattened_embeddings = embeddings.reshape(embeddings.shape[0], -1).transpose()

        logger.info('Clustering embeddings...')
        # perform clustering and reshape in order to get the segmentation volume
        start = time.time()
        clusters = self.clustering.fit_predict(flattened_embeddings).reshape(output_shape)
        logger.info(
            f'Number of clusters found by {self.clustering}: {np.max(clusters)}. Duration: {time.time() - start} sec.')
        return clusters

    def _merge_segmentation(self, segmentation, index, output_segmentation, visited_voxels_array):
        """
        Given the `segmentation` patch, its `index` in the `output_segmentation` array and the array visited voxels
        merge the segmented patch (`segmentation`) into the `output_segmentation`

        Args:
            segmentation (ndarray): segmented patch
            index (tuple): position of the patch inside `output_segmentation` volume
            output_segmentation (ndarray): current state of the output segmentation
            visited_voxels_array (ndarray): array of voxels visited so far (same size as `output_segmentation`); visited
                voxels will be marked by a number greater than 0
        """
        index = tuple(index)
        # get new unassigned label
        max_label = np.max(output_segmentation) + 1
        # make sure there are no clashes between current segmentation patch and the output_segmentation
        # but keep the noise label
        noise_mask = segmentation == self.noise_label
        segmentation += int(max_label)
        segmentation[noise_mask] = self.noise_label
        # get the overlap mask in the current patch
        overlap_mask = visited_voxels_array[index] > 0
        # get the new labels inside the overlap_mask
        new_labels = np.unique(segmentation[overlap_mask])
        merged_labels = self._merge_labels(output_segmentation[index], new_labels, segmentation)
        # relabel new segmentation with the merged labels
        for current_label, new_label in merged_labels:
            segmentation[segmentation == new_label] = current_label
        # update the output_segmentation
        output_segmentation[index] = segmentation
        # visit the patch
        visited_voxels_array[index] += 1

    def _merge_labels(self, current_segmentation, new_labels, new_segmentation):
        def _most_frequent_label(labels):
            unique, counts = np.unique(labels, return_counts=True)
            ind = np.argmax(counts)
            return unique[ind]

        result = []
        # iterate over new_labels and merge regions if the IoU exceeds a given threshold
        for new_label in new_labels:
            # skip 'noise' label assigned by hdbscan
            if new_label == self.noise_label:
                continue
            new_label_mask = new_segmentation == new_label
            # get only the most frequent overlapping label
            most_frequent_label = _most_frequent_label(current_segmentation[new_label_mask])
            # skip 'noise' label
            if most_frequent_label == self.noise_label:
                continue
            current_label_mask = current_segmentation == most_frequent_label
            # compute Jaccard index
            iou = np.bitwise_and(new_label_mask, current_label_mask).sum() / np.bitwise_or(new_label_mask,
                                                                                           current_label_mask).sum()
            if iou > self.iou_threshold:
                # merge labels
                result.append((most_frequent_label, new_label))

        return result

    def _get_clustering(self, clustering_alg, kwargs):
        logger.info(f'Using {clustering_alg} for clustering')

        if clustering_alg == 'hdbscan':
            min_cluster_size = kwargs.get('min_cluster_size', 50)
            min_samples = kwargs.get('min_samples', None),
            metric = kwargs.get('metric', 'euclidean')
            cluster_selection_method = kwargs.get('cluster_selection_method', 'eom')

            logger.info(f'HDBSCAN params: min_cluster_size: {min_cluster_size}, min_samples: {min_samples}')
            return hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples, metric=metric,
                                   cluster_selection_method=cluster_selection_method)
        else:
            bandwidth = kwargs['bandwidth']
            logger.info(f'MeanShift params: bandwidth: {bandwidth}, bin_seeding: True')
            # use fast MeanShift with bin seeding
            return MeanShift(bandwidth=bandwidth, bin_seeding=True)
