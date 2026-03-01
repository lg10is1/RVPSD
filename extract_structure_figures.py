
import os
from pymol import cmd

# start PyMOL
cmd.reinitialize()


pdb_folder = "$path_to/pdb_files"
output_folder = "$path_to/all_pdb_figure_png"


pdb_files = [f for f in os.listdir(pdb_folder) if f.endswith('.pdb')]

cmd.set_color("protein_blue", [0.0, 0.0, 1.0]) 


for pdb_file in pdb_files:
    
    cmd.load(os.path.join(pdb_folder, pdb_file))
    
    
    cmd.color("protein_blue", "all")
    
    cmd.set("ray_opaque_background", 0)
    
    
    output_filename = os.path.splitext(pdb_file)[0] + ".png"  
    cmd.png(os.path.join(output_folder, output_filename), width=330, height=330, dpi=300, ray=1)  
    
    
    cmd.delete('all')

# quit PyMOL
cmd.quit()