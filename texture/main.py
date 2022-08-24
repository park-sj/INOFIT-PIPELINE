# -*- coding: utf-8 -*-
"""
Created on Thu May 20 16:29:54 2021

@author: user
"""

import sys
import os
import bpy
import matplotlib.pyplot as plt

def execute_quiet(func, *args, **kwargs):
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, "w")
    result = func(*args, **kwargs)
    sys.stdout = old_stdout
    return result


def bake_texture(scanner_loc, innofit_loc, save_dir):
    # First delete all the default objects - cube, camera, light
    bpy.ops.object.select_all(action='SELECT')
    execute_quiet(bpy.ops.object.delete)

    # 입력 받은 해상도에 맞춰서 저장해줘야해서 png 열어서 해상도 확인
    # png 파일 확인
    files = os.listdir(os.path.dirname(scanner_loc))
    for file in files:
        if os.path.splitext(file)[1] == '.png':
            png_path = os.path.join(os.path.dirname(scanner_loc), file)
#    png_path = scanner_loc[:-4] + '_0.png'
    resolution = plt.imread(png_path, format='png').shape

    ## Import innofit mesh and prepare to bake
    # 1. import
    imported_object = execute_quiet(bpy.ops.import_scene.obj, filepath=innofit_loc, split_mode='OFF')
    assert 'FINISHED' in imported_object, 'Failed to load innofit mesh!'
    innofit_obj = bpy.data.objects[0]    
    # 2. unwrap
    bpy.context.view_layer.objects.active = innofit_obj
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.uv.unwrap()
    bpy.ops.object.mode_set(mode='OBJECT')
    print('Innofit mesh is unwrapped')
    # 3. make and connect nodes
    for mat in innofit_obj.data.materials:
        print(f'mat: {mat}')
        mat.use_nodes = True
        nodes = mat.node_tree.nodes
        links = mat.node_tree.links
        ps = nodes.get('Principled BSDF')
        # 3-a. Add image texture node and make new image
        texture_node = nodes.new('ShaderNodeTexImage')
        texture_node.name = 'Image Texture'
        texture_node.select = True
        nodes.active = texture_node
        image_name = innofit_obj.name + '_BakedTexture'
        img = bpy.data.images.new(image_name, resolution[0], resolution[1])
        texture_node.image = img
        links.new(texture_node.outputs['Color'], ps.inputs['Base Color'])
        # 3-b. Add Mapping node
        mapping_node = nodes.new('ShaderNodeMapping')
        mapping_node.select = True
        links.new(mapping_node.outputs['Vector'], texture_node.inputs['Vector'])
        # 3-c. Add UV Map node
        uvmap_node = nodes.new('ShaderNodeUVMap')
        uvmap_node.select = True
        links.new(uvmap_node.outputs['UV'], mapping_node.inputs['Vector'])

    # import scanner mesh
    imported_object = execute_quiet(bpy.ops.import_scene.obj, filepath=scanner_loc)
    assert 'FINISHED' in imported_object
    for obj in bpy.data.objects:
        if obj != innofit_obj:
            scanner_object = obj
    # scanner_object = bpy.data.objects[1]

    # bake
    print("Selected objects: ", bpy.context.selected_objects)
    bpy.context.scene.render.engine = 'CYCLES'
    
    scanner_object.select_set(True)
    '''
    bpy.context.scene.cycles.bake_type = 'DIFFUSE'
    bpy.context.scene.render.bake.use_pass_direct = False
    bpy.context.scene.render.bake.use_pass_indirect = False
    bpy.context.scene.render.bake.use_selected_to_active = True
    bpy.context.scene.render.bake.max_ray_distance = 0.02
    bpy.ops.object.bake_image()
    '''
    bpy.ops.object.bake(type="DIFFUSE", pass_filter={'COLOR'}, use_selected_to_active=True, target='IMAGE_TEXTURES', max_ray_distance=1, cage_extrusion=0.5)

    # save results
    img.save_render(filepath=os.path.join(save_dir, 'result_trimmed_baked.png'))
    
    bpy.ops.object.select_all(action='DESELECT')
    scanner_object.select_set(True)
    bpy.ops.object.delete()

    innofit_obj.select_set(True)

    output_loc = os.path.join(save_dir, 'result_trimmed_baked.obj')
    execute_quiet(bpy.ops.export_scene.obj, filepath=output_loc, use_uvs=True, use_normals=False, use_materials=True, keep_vertex_order=True, use_selection=True)
    
    with open(os.path.join(save_dir, 'result_trimmed_baked.mtl'), 'r') as f:
        filedata = f.read()
    filedata = filedata.replace('map_Kd .', 'map_Kd result_trimmed_baked.png')

    with open(os.path.join(save_dir, 'result_trimmed_baked.mtl'), 'w') as f:
        f.write(filedata)
    
if __name__ == "__main__":
    if len(sys.argv) != 4:
        scanner_loc = "/home/shkim/Libraries/texture_transfer/data/sb.obj"
        innofit_loc = "/home/shkim/Libraries/texture_transfer/data/sb_innofit.obj"
        save_dir = '/home/shkim/Libraries/texture_transfer/data/'
    else:
        scanner_loc = sys.argv[1]
        innofit_loc = sys.argv[2]
        save_dir = sys.argv[3]
    assert os.path.isfile(scanner_loc), f'{scanner_loc} file does not exist!'
    assert os.path.isfile(innofit_loc), f'{innofit_loc} file does not exist!'

    bake_texture(scanner_loc, innofit_loc, save_dir)
