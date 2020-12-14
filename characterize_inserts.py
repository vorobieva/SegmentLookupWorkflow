import pyrosetta as py
import glob
import re
from shutil import copyfile
import sys,copy
import os
import numpy as np
from math import sin,cos,pi,asin
from collections import defaultdict
import pandas as pd

metal_coord = [-7.77678571,  6.31519048,  9.5682619 ]

py.init()

def read_PDBInfo(pdb,label):
	label_res = []
	with open(pdb, 'r') as in_pdb:
		for line in in_pdb:
			if "REMARK PDBinfo-LABEL:" in line and line.split()[3] == label:
				label_res.append(int(line.split()[2]))
	return label_res

def inproduct(v1,v2):
	sum=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
	return(sum)

def cross(v1,v2):
	v=[v1[1]*v2[2]-v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0]]
	return(v)

def angle(v1,v2):
	inprod=0
	for i in range(3):
		inprod+=v1[i]*v2[i]
	dist1=(v1[0]**2+v1[1]**2+v1[2]**2)**0.5
	dist2=(v2[0]**2+v2[1]**2+v2[2]**2)**0.5
	sinval=inprod/(dist1*dist2)
	if sinval >1.0 or sinval <-1.0:
		angle=0
	else:
		angle=asin(sinval)-pi/2.0
	return(angle)

def pdb2crd(pdbfile,atmname):
	pdbcont=file(pdbfile)
	crd = {}
	for line in pdbcont:
		if not line.startswith("ATOM"):
			continue
		resno = int(line[22:26])
		atmtype = line[12:16].strip()
		if atmtype == atmname:
			crd[resno] = [float(line[30+i*8:38+i*8]) for i in range[3]]
	return(crd)

def pdb2contacts(pdbfile,ref_coord,dist,scaffold_contacts=None,loop_params=None):
	pdbcont=open(pdbfile,"r")
	contacts = []
	heavy_atms = ["CA","C","N","O","CB"]
	for line in pdbcont:
		if not line.startswith("ATOM"):
			continue
		coord = [line.split()[6],line.split()[7],line.split()[8]]
		v = [float(coord[k]) - ref_coord[k] for k in range(3)]
		norm = inproduct(v,v)**0.5
		if norm <= dist and line.split()[2] in heavy_atms:
			atm = line.split()[2] + "_" + line.split()[5] + "_" + line.split()[4]
			if scaffold_contacts==None and loop_params==None:
				contacts.append(atm)
			elif scaffold_contacts!=None and loop_params!=None:
				reindexed_ref_contacts = []
				for pos in scaffold_contacts:
					vals = pos.split("_")
					if int(vals[1]) <= loop_params[1]+1:
						reindexed_ref_contacts.append(pos)
					elif int(vals[1]) > loop_params[1]+1 and int(vals[1]) < (loop_params[2]+loop_params[0]-1):
						continue
					elif int(vals[1]) >= (loop_params[2]+loop_params[0]-1):
						new_val =  int(vals[1]) + loop_params[0]
						reindexed_ref_contacts.append("_".join([vals[0],str(new_val),vals[2]]))
					else:
						print("Something is wrong with the indexing!")
				if not atm in reindexed_ref_contacts:
					contacts.append(atm)
	pdbcont.close()
	return(contacts)

def CBcontacts(pdbfile,contacts,metal_coord):
	pdbcont=open(pdbfile,"r")
	contacts_coord = defaultdict(dict)
	contact_res = []
	contact_res_CB = []
	CB_contacts = []
	heavy_atms = ["CA","C","N","O","CB"]
	for contact in contacts:
		if contact.split("_")[1] not in [sl[0] for sl in contact_res]:
			contact_res.append([contact.split("_")[1],contact.split("_")[2] ])
		if contact.split("_")[0] == "CB":
			contact_res_CB.append([contact.split("_")[1],contact.split("_")[2]])
	for line in pdbcont:
		if not line.startswith("ATOM"):
			continue
		resno =  line.split()[5]
		chain_id = line.split()[4]
		atmtype =  line.split()[2]
		if [resno, chain_id] in contact_res_CB and atmtype in heavy_atms:
			res_id = resno + chain_id
			contacts_coord[res_id][atmtype] =  [float(line.split()[6]),float(line.split()[7]),float(line.split()[8])]
	for key in contacts_coord:
		v_ref = [contacts_coord[key]["CA"][l] - metal_coord[l] for l in range(3)]
		norm_ref = inproduct(v_ref,v_ref)**0.5
		v_ref = [a/norm_ref for a in v_ref]
		v = [contacts_coord[key]["CB"][l] - contacts_coord[key]["CA"][l] for l in range(3)]
		norm = inproduct(v,v)**0.5
		v = [a/norm for a in v]
		if inproduct(v_ref,v) <= -0.4:
			CB_contacts.append(key)
	pdbcont.close()
	return(CB_contacts)

#rms = py.rosetta.protocols.analysis.simple_metrics.RMSDMetric()
abego = py.rosetta.core.sequence.ABEGOManager()
DSSP = py.rosetta.protocols.moves.DsspMover()
scorefxn = py.create_score_function("fldsgn_cen.wts")
scorefxn.set_weight(py.rosetta.core.scoring.hbond_lr_bb, 1.0) 
scorefxn.set_weight(py.rosetta.core.scoring.hbond_sr_bb, 1.0) 
cen = py.rosetta.protocols.simple_moves.SwitchResidueTypeSetMover( 'centroid' )

cols=["template","scaffold","loop number","loop secondary structure","loop abego string","loop length","loop start position","loop end position","protein length","fraction of L residues","fraction of H residues","fraction of E residues","longest stretch of L residues","long range backbone hydrogen bond energy","short range backbone hydrogen bond energy","backbone clashes (vdw)","RMSD to scaffold PDB", "Atoms at 8A from metal","Atoms at 9A from metal","Atoms at 10A from metal","Atoms at 11A from metal","Atoms at 12A from metal","Number atoms at 8A","Number atoms at 9A","Number atoms at 10A","Number atoms at 11A","Number atoms at 12A","CB contacts at 8A","CB contacts at 9A","CB contacts at 10A","CB contacts at 11A","CB contacts at 12A","Number CB contacts at 8A","Number of CB contacts at 9A","Number of CB contacts at 10A","Number of CB contacts at 11A","Number of CB contacts at 12A"]
df = pd.DataFrame(columns=cols)

for task_file in glob.glob("scaffolds_A/task_list"):
	folder = os.path.dirname(task_file)
	scaffold = folder.split("_")[1]
	ref_pdb = folder+ "/templates/ref.pdb"
	pose_ref = py.pose_from_pdb(ref_pdb)
	cen.apply(pose_ref)
	scorefxn(pose_ref)
	hb_lr_ref = pose_ref.energies().total_energies()[py.rosetta.core.scoring.hbond_lr_bb]
	hb_sr_ref = pose_ref.energies().total_energies()[py.rosetta.core.scoring.hbond_sr_bb]
	l_ref = pose_ref.total_residue()
	pdb_ref_dimer = folder+ "/templates/ref_dimer.pdb"
	contacts_ref = {8:[], 9:[], 10:[], 11:[], 12:[]}
	for d in range(8,13):
		contacts = pdb2contacts(pdb_ref_dimer,metal_coord,d,scaffold_contacts=None,loop_params=None)
		contacts_ref[d] = contacts		
	with open(task_file) as task:
		for line in task:
			pdb_file = folder + "/" + line.strip()
			pdb_name= os.path.splitext(line.strip())[0]
			pdb_dimer = folder + "/" + pdb_name + "_0001.pdb"
			template = pdb_name.split("_")[1]
			loop_no = pdb_name.split("_")[-1]
			pose = py.pose_from_pdb(pdb_file)
			cen.apply(pose)
			scorefxn(pose)
			rms = py.rosetta.core.scoring.CA_rmsd(pose,pose_ref)
			DSSP.apply(pose)
			ss = pose.secstruct()
			hb_lr = pose.energies().total_energies()[py.rosetta.core.scoring.hbond_lr_bb]
			hb_sr = pose.energies().total_energies()[py.rosetta.core.scoring.hbond_sr_bb]
			lr = hb_lr - hb_lr_ref
			sr = hb_sr - hb_sr_ref
			vdw = pose.energies().total_energies()[py.rosetta.core.scoring.vdw]
			segment = read_PDBInfo(pdb_file,"segment_lookup")
			loop_len = len(segment)
			loop_start = min(segment)
			loop_end = max(segment)
			b = abego.get_symbols(pose)
			a = abego.get_abego_string(b)
			l = pose.total_residue()
			added_len = l - l_ref
			loop_ss = ""
			loop_abego = ""
			for resi in segment: 
				loop_ss += ss[resi-1]
				loop_abego += a[resi-1]
			list_of_L_segments = re.findall(r"[L]+", loop_ss)
			list_of_H_segments = re.findall(r"[H]+", loop_ss)
			list_of_E_segments = re.findall(r"[H]+", loop_ss)
			L_content = len("".join(list_of_L_segments))/loop_len
			H_content = len("".join(list_of_H_segments))/loop_len
			E_content = len("".join(list_of_E_segments))/loop_len
			longest_L_stretch = max(list_of_L_segments, key=len)
			contacts_dict = {8:[], 9:[], 10:[], 11:[], 12:[]}
			CB_dict = {8:[], 9:[], 10:[], 11:[], 12:[]}
			for d in range(8,13):
				contacts = pdb2contacts(pdb_dimer,metal_coord,d,scaffold_contacts=contacts_ref[d],loop_params=[added_len,loop_start,loop_end])
				contacts_dict[d] = contacts
				CB_pointing_to_ligand = CBcontacts(pdb_dimer,contacts,metal_coord)
				CB_dict[d] = CB_pointing_to_ligand

			data = pd.Series([template,scaffold,loop_no,loop_ss,loop_abego,loop_len,loop_start,loop_end,l,L_content,H_content,E_content,longest_L_stretch,lr,sr,vdw,rms,",".join(contacts_dict[8]),",".join(contacts_dict[9]),",".join(contacts_dict[10]),",".join(contacts_dict[11]),",".join(contacts_dict[12]),len(contacts_dict[8]),len(contacts_dict[9]),len(contacts_dict[10]),len(contacts_dict[11]),len(contacts_dict[12]),",".join(CB_dict[8]),",".join(CB_dict[9]),",".join(CB_dict[10]),",".join(CB_dict[11]),",".join(CB_dict[12]),len(CB_dict[8]),len(CB_dict[9]),len(CB_dict[10]),len(CB_dict[11]),len(CB_dict[12])], index=cols )
			df = df.append(data, ignore_index=True)

df.to_pickle("segment_data_A.pkl")

#	if len(longest_L_stretch) < 13:
#		start_loop = min(segment)
#		stop_loop = max(segment)
#		stop_scaffold1 = start_loop - 1
#		start_scaffold2 = stop_loop + 1
#		stop_scaffold2 = len(ss)
#		copyfile(pdb_file, "hybridize/"+pdb_file)
#		with open("hybridize/task_list",'a') as slurm_jobs:
#			slurm_jobs.write("%s %d %d %d %d %d\n" %(pdb_file, stop_scaffold1, start_loop, stop_loop, start_scaffold2, stop_scaffold2))
