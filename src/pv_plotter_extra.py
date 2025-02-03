#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 10:11:31 2024

@author: hector.ortiz
"""

import numpy as np
import pandas as pd
import pyvista as pv
import re
from forces_eq import *

def save_meshes(robot,anim_data):
    robot.OUTPUT_FLAG = True
    print('Saving Meshes')
    #for i in range():
    N = anim_data.shape[0]
    for i in range(N):    
        #os.system('cls' if os.name == 'nt' else 'clear')
        if i==N:
            print(f'saving mesh: {i}/{N}') 
        else:
            print(f'saving mesh: {i}/{N}', end='\r') 
        
        #r1,r2,r3 = y[i].copy()
        R = anim_data[i,6:]
        X = anim_data[i,:6]
        #robot.forward_kinematics(r1,r2,r3)
        robot.eqsystem_forward(X, *R)
        
        robot.save_mesh(i)
    robot.OUTPUT_FLAG = False
    print('All meshed saved')

def update_indices(index_array, data_array):
    # Create a mapping of old indices to new indices
    unique_indices = np.unique(index_array)
    index_map = {old_idx: new_idx for new_idx, old_idx in enumerate(unique_indices)}
    
    # Apply the mapping to the first three columns of the data_array
    updated_data = data_array.copy()
    
    # Replace old indices with new indices
    for col in range(3):  # First three columns
        updated_data[:, col] = np.vectorize(index_map.get)(data_array[:, col])
    
    return updated_data

def save_droplet_stl():
    vertex_l, facets_l = get_list_vert_facet_files()
    print(vertex_l)
    print(facets_l)

    # Extract the part before '.txt'
    for filename in vertex_l:
        base = filename.split(".txt")[0]
        # Extract the last two parts split by "_"
        parts = base.split("_")[-2:]  # Take only the last two parts
        # Extract i and j
        i, j = map(int, parts)
        # Get mesh
        mesh,c,a = mesh_droplet(j,i)["mesh0"]
        mesh.save(f'../data/mesh/mesh_{i}_{j}.stl')

def get_list_vert_facet_files():
        # Folder containing the files
    folder_path = '../data/mesh'

    # Lists to hold filenames
    vertex_files = []
    face_files = []

    # Dictionary to group by (i, j)
    grouped_files = {}

    # Iterate through files in the folder
    for filename in os.listdir(folder_path):
        if filename.startswith("vertex_") or filename.startswith("facet_"):
            parts = filename.split("_")
            if len(parts) == 3:
                prefix, i, j_with_ext = parts
                j = os.path.splitext(j_with_ext)[0]  # Remove extension from j
                key = (i, j)
                
                # Initialize the key in the dictionary if it doesn't exist
                if key not in grouped_files:
                    grouped_files[key] = {"vertex": None, "facet": None}
                
                # Assign filenames to the correct group
                if prefix == "vertex":
                    grouped_files[key]["vertex"] = filename
                elif prefix == "facet":
                    grouped_files[key]["facet"] = filename

    # Separate the grouped files into respective lists
    for key, files in grouped_files.items():
        if files["vertex"]:
            vertex_files.append(files["vertex"])
        if files["facet"]:
            face_files.append(files["facet"])
    return vertex_files, face_files


def mesh_droplet(i,j,w=[0,0,0], save_stl=False):
    
    file_path = f'../data/mesh/vertex_{j}_{i}.txt'
    data = pd.read_csv(file_path, sep=' ', header=None).to_numpy()
    
    #verts = data[:,1:]*1e6
    verts = data[:,:].copy()
    index_array = verts[:,0].astype(int)
    verts[:,1:] = verts[:,1:]*1e6
    verts[:,1] = verts[:,1] + w[0]
    verts[:,2] = verts[:,2] + w[1]
    
    points = verts[:,1:].copy()
    
    data = np.loadtxt(f'../data/mesh/facet_{j}_{i}.txt')
    facets = data.copy().astype(int)
    face_array = facets[:,1:]
    mapping = {original: new for new, original in enumerate(index_array)}
    
    faces_orig = face_array.copy()  # Keep the original reference array intact
    faces_orig[:, :3] = np.vectorize(mapping.get)(faces_orig[:, :3])

    val = ['blue','red','black', 'white']
    
    # Get unique values in the last column
    unique_values = np.unique(faces_orig[:, -1])
    
    # Separate rows by the value in the last column
    #ft = [faces_orig[faces_orig[:, -1] == value,:-1] for value in unique_values]
    ft =  [faces_orig[faces_orig[:, 3] == value,:-1] for value in unique_values]

    f1 = ft[0]
    f2 = ft[1]
    f3 = ft[2]
    f = [f1,f2,f3]
    
    faces2 = np.hstack((np.full((f[0].shape[0], 1), 3),f[0]))
    faces2 = faces2[:3,:]
    
    #mesh0 = pv.PolyData(points,  np.hstack((np.full((f[1].shape[0], 1), 3),f[1])))
    mesh0 = pv.make_tri_mesh(points,  f[0])
    mesh1 = pv.make_tri_mesh(points,  f[1])
    mesh2 = pv.make_tri_mesh(points,  f[2])
    color0 = "cornflowerblue"
    alpha0 = 0.7
    color1 = "red"
    alpha1 = 1
    
    # Create a dictionary to hold meshes and their colors
    mesh_d = {
               "mesh0": (mesh0, color0, alpha0),
             #  "mesh1": (mesh1, color1, alpha1)
              }
    if save_stl:
        base_path = f'../data/mesh/mesh_{i}_{j}'
        file_path = f'{base_path}.stl'
        counter = 1
        
        # Check if file exists and modify the name until it's unique
        while os.path.exists(file_path):
            file_path = f'{base_path}_{counter}.stl'
            counter += 1
        
        mesh0.save(file_path)
    return mesh_d

def create_mesh_stl(file_path, P,u_r, scale=1e3, save_stl=False): 
    mesh = pv.read(file_path).scale(scale)
    R = R_from_vec(u_r)
    M_44 = np.zeros((4,4))
    M_44[:-1,-1] = P#/1000
    M_44[:-1,:-1] = R
    M_44[-1,-1] = 1
    return mesh.transform(M_44)

def plot_robot(plotter,i,w,P,u_r,file_path, save_stl=False):
    for j in range(1,4):
        mesh_d = mesh_droplet(i, j, w[j-1], save_stl)
        for name, (mesh, color, opacity) in mesh_d.items():
            _=plotter.add_mesh(mesh, color=color, opacity=opacity, label=name,
                               split_sharp_edges=True,
                               diffuse=0.6, specular=0.4, ambient=0.5)
    # Plot tripod

    mesh = create_mesh_stl(file_path, P, u_r)
    if save_stl:
        base_path = f'../data/mesh/platform_{i}_{j}'
        file_path = f'{base_path}.stl'
        counter = 1
        
        # Check if file exists and modify the name until it's unique
        while os.path.exists(file_path):
            file_path = f'{base_path}_{counter}.stl'
            counter += 1
        
        mesh.save(file_path)
    
    # Plot table
    vertices = np.array([[0, 0, 0], [2, 0, 0], [2, 2, 0], [0, 2, 0]]).astype(float)*1000.0*1.5
    vertices[:,:2] = vertices[:,:2] - 1000*1.5
    faces = [4,0,1,2,3]
    mesh_table = pv.PolyData(vertices, faces)
    _= plotter.add_mesh(mesh_table, color='white',
                    #    pbr=True,
                        specular=0.5, diffuse=0.5,ambient=0.5,specular_power=30)
    
    _= plotter.add_mesh(mesh, color="red", show_edges=False, split_sharp_edges=True,
                     #   pbr=True,
                        specular=0.5, diffuse=0.5,ambient=0.5,specular_power=30)


if __name__ == "__main__":
    save_droplet_stl()
