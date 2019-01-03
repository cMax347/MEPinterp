import numpy as np

#   this script
#           READS the 'kpts' file generated with kptsgen.pl (i.e. FLEUR kpt-list format)
#       &
#           WRITES the 'seed_name_geninterp.kpt' file in current directory
#
#
#
#     note: all files should be present within the source dir




#read header (1.line)
with open('./kpts','r') as kpt_file:
    first_line  =   kpt_file.readline()

#get 1. line parameters
paras       =   np.fromstring(first_line, dtype=float, count=2, sep=' ')
glob_scal   =   paras[1]
nkpts_para  =   int( paras[0] )




#get the kpt lst (read body)
kpt_lst     =   np.loadtxt('./kpts',skiprows=1, usecols=(0,1,2,3) )
nKpts       =   kpt_lst.shape[0]




#   print para info
if nKpts != nkpts_para:
    print("WARNING expectec"+str(nkpts_para)+' kpts, but found'+str(nKpts))
else:
    print('found '+str(nKpts)+' kpts as expected')

print('global scaling factor: '+'{0:.8e}'.format(glob_scal) )




#   WRITE OUTPUT

seed_name   =   'lead'

outFile =   './'+seed_name+'_geniterp.kpt'
with open(outFile ,'w') as out:
    #   WRITE HEADER
    out.write('#random comment\n')
    out.write('frac\n')
    out.write(str(nKpts)+'\n')

    #   WRITE KPTS
    for line, kpt in enumerate(kpt_lst):
        out.write(str(line+1)+' '+str(glob_scal*kpt[0]*kpt[3])+' '+str(glob_scal*kpt[1]*kpt[3])+' '+str(glob_scal*kpt[2]*kpt[3])+'\n')
close(out)
print("wrote "+outFile+" , by!")

