import random as rd
import os
import os.path

file_name = "cc"
dir_ = os.path.dirname(os.path.abspath(file_name)) + "\\"
f_cc_name = dir_ + file_name + ".txt"
f_cc= open(f_cc_name,"w+")
lead_str = "    82    -11.5 "
poly_str = "    112   -0.93 "
Imp_str = "  Imp:p,n=1 \n"
#SHELLS:
#Shell 1: 1-240, 15*16=240
#Shell 2:  1001-1345, 15*23=345
#Shell 3:  2001-2465, 15*31=465
#Shell 4:  3001-3570, 15*38=570

#Configuration to Generate:
#configuration = 'Checkerboard'
configuration = 'Random'

if (configuration == 'Checkerboard'):
    flip=False
    for cell in range(240):
        if (cell%16==0):
            flip = not flip
        if flip:
            if (cell % 2) == 0:
                cellcard = str(cell+1)+lead_str+str(-cell-1) + "\n" #+Imp_str
            else:
                cellcard = str(cell+1)+poly_str+str(-cell-1) + "\n" #+Imp_str
            f_cc.write(cellcard)
        else:
            if (cell % 2) == 0:
                cellcard = str(cell+1)+poly_str+str(-cell-1) + "\n"  #+Imp_str
            else:
                cellcard = str(cell+1)+lead_str+str(-cell-1) + "\n"  #+Imp_str
            f_cc.write(cellcard)
    for cell in range(1000,1345):
        if (cell % 2) == 0:
            cellcard = str(cell+1)+poly_str+str(-cell-1) + "\n"  #+Imp_str
        else:
            cellcard = str(cell+1)+lead_str+str(-cell-1) + "\n"  #+Imp_str
        f_cc.write(cellcard)    
    for cell in range(2000,2465):
        if (cell % 2) == 0:
            cellcard = str(cell+1)+lead_str+str(-cell-1) + "\n"  #+Imp_str
        else:
            cellcard = str(cell+1)+poly_str+str(-cell-1) + "\n"  #+Imp_str
        f_cc.write(cellcard)  
    for cell in range(3000,3570):
        if (cell%38==0):
            flip = not flip
        if flip:
            if (cell % 2) == 0:
                cellcard = str(cell+1)+poly_str+str(-cell-1) + "\n"  #+Imp_str
            else:
                cellcard = str(cell+1)+lead_str+str(-cell-1) + "\n"  #+Imp_str
            f_cc.write(cellcard)
        else:
            if (cell % 2) == 0:
                cellcard = str(cell+1)+lead_str+str(-cell-1) + "\n"  #+Imp_str
            else:
                cellcard = str(cell+1)+poly_str+str(-cell-1) + "\n"  #+Imp_str
            f_cc.write(cellcard)
elif (configuration == 'Random'):
    for cell in range(240):
        randnum = rd.uniform(0,1)
        if (randnum >= 0.5) == 0:
            cellcard = str(cell+1)+lead_str+str(-cell-1) + "\n"  #+Imp_str
        else:
            cellcard = str(cell+1)+poly_str+str(-cell-1) + "\n"  #+Imp_str
        f_cc.write(cellcard)
    for cell in range(1000,1345):
        randnum = rd.uniform(0,1)
        if (randnum >= 0.5) == 0:
            cellcard = str(cell+1)+lead_str+str(-cell-1) + "\n"  #+Imp_str
        else:
            cellcard = str(cell+1)+poly_str+str(-cell-1) + "\n"  #+Imp_str
        f_cc.write(cellcard)    
    for cell in range(2000,2465):
        randnum = rd.uniform(0,1)
        if (randnum >= 0.5) == 0:
            cellcard = str(cell+1)+lead_str+str(-cell-1) + "\n"  #+Imp_str
        else:
            cellcard = str(cell+1)+poly_str+str(-cell-1) + "\n"  #+Imp_str
        f_cc.write(cellcard)  
    for cell in range(3000,3570):
        randnum = rd.uniform(0,1)
        if (randnum >= 0.5) == 0:
            cellcard = str(cell+1)+lead_str+str(-cell-1) + "\n"  #+Imp_str
        else:
            cellcard = str(cell+1)+poly_str+str(-cell-1) + "\n"  #+Imp_str
        f_cc.write(cellcard)

f_cc.close()