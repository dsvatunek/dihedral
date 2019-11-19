# select atoms based on numbers starting with 1!
selected = [8,1,4,7]

selected = [int(x)-1 for x in selected]

class structures:
	pass
    
import glob
import numpy as np

periodic_table = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
    "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]
	
def isInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False	
        
def getDihedral(p0, p1, p2, p3):
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))
    
    
def structures_from_xyz(file):
	structures.atoms = []
	structures.xyz = []
	structures.title = []
	file_object = open(file, 'r')
	input = (line for line in file_object) #make generator
	#search for number of atoms
	for line in input:
		if isInt(line.strip()):
			n_atoms=int(line)
			break	
	else: #exits if no line with number of atoms was found
		sys.exit('Error:\t\tNo xyz coordinates found in file: ' + file)
	#skip one line
	structures.title.append(next(input).strip())
	# now there should be n_atoms lines of coordinates
	for i in range(n_atoms):
		l=next(input).split()
		if l[0] in periodic_table:
			structures.atoms.append(l[0]) #get atom symbol and append to atom list
		else:
			sys.exit('Error:\t\tsomething is wrong with the first structure in file: '+file)
		coords=[float(x) for x in l[1:]] #convert line to list of floats
		coords=np.array([coords]) #create array with coords
		try: #try append, doesn't work if XYZ doesn't exist yet
			XYZ=np.concatenate((XYZ,coords), axis=0)
		except NameError:
			XYZ=coords
	structures.xyz.append(XYZ) #append first structure to structures list
	del XYZ #get rid of that for the next structure
	#now search for more structures
	for line in input:
		#start extracting if atom number line is found
		try:
			if int(line.strip()) == n_atoms:
				#read one line to skip title
				structures.title.append(next(input).strip())
				# now there should be n_atoms lines of coordinates
				for i in range(n_atoms):
					l=next(input).split()
					coords=[float(x) for x in l[1:]]
					coords=np.array([coords])
					try: #try append, doesn't work if XYZ doesn't exist yet
						XYZ=np.concatenate((XYZ,coords), axis=0)
					except NameError:
						XYZ=coords
				structures.xyz.append(XYZ)
				del XYZ
		except ValueError:
			pass			
	return structures
    
    
def main():
    filelist=sorted(glob.glob("*.xyz"))
    for file in filelist:
        structures = structures_from_xyz(file)
        output = open(file+'_dihedrals.csv',"w")
        for i in range(len(structures.xyz)):
            output.write(str(getDihedral(structures.xyz[i][selected[0]],structures.xyz[i][selected[1]],structures.xyz[i][selected[2]],structures.xyz[i][selected[3]]))+"\n")
        output.close()
            
if __name__ == "__main__":
    main()
