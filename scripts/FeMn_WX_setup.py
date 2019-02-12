import os
import datetime
import numpy as np
import shutil as sh
from FeMn_latt_OPT import FeMn_latt_OPT
from tb_input_writer import write_mepInterp_input



# Pauli matrix
sigma_x = np.array([[0 ,  1], [1 , 0]])
sigma_y = np.array([[0 ,-1j], [1j, 0]])
sigma_z = np.array([[1 ,  0], [0 ,-1]])


# SPIN configs (radiants)# (theta,phi) are given in radians, which correspond to (beta,alpha) in Fleur.
theta_3q = [ 0.0000000000,  1.9106332362,  1.9106332362, 1.9106332362]
phi_3q   = [ 0.0000000000, -2.0943951024,  2.0943951024, 0.0000000000]

theta_2q = [ 0.6154797087,  2.5261129449,  1.5707963268,  1.5707963268]
phi_2q   = [ 1.0471975512, -2.0943951024,  2.6179938780, -0.5235987756]

theta_1q = [ 0.9553166181,  0.9553166181,  2.1862760355,  2.1862760355]
phi_1q   = [-2.0943951024, -2.0943951024,  1.0471975512,  1.0471975512]

theta_0q = [ 0.0000000000,  1.5707963268,  1.5707963268,  1.5707963268] # For test
phi_0q   = [ 0.0000000000, -2.0943951024,  2.0943951024,  0.0000000000]


#init latt
latt0 = np.array([
            [  5.60092859 , 0.00000000,  3.96045459],   # original lattice
            [ -2.80046430 , 4.85054644,  3.96045459],  # in unit of Bohr
            [ -2.80046430 ,-4.85054644,  3.96045459]])
vol = np.dot(   latt0[0,:], np.cross(latt0[1,:],latt0[2,:]))


def print_v(msg,v):
    if v:
        print(msg)


def dir_setup(base_dir, sub_dir,verbose=False):
    #
    #   force correct dir string format
    if not (base_dir[-1] == '/'):
        base_dir    =   base_dir+'/'
    if not (sub_dir[-1] == '/'):
        sub_dir    =   sub_dir+'/'
    #
    #
    #   make sure the base folder exists
    if not os.path.isdir(base_dir):
        os.mkdir(base_dir)
    #
    #   make sure root dir is clean
    root_dir     =   base_dir+sub_dir
    if os.path.isdir(root_dir):
        sh.rmtree(root_dir)
        print_v('[FeMn_WX_setup]: removed old target folder: '+root_dir, verbose)
    os.mkdir(root_dir)
    #
    #   create the w90 file container dir
    w90_dir      =   root_dir+'w90files/'
    os.mkdir(w90_dir)
    #
    return root_dir, w90_dir


def write_3q_HR(    base_dir,sub_dir, seed_name, t,strain, J_ex,Jz,tso, spin_order,
                    mp_grid, kubo_tol, valence_bands,
                    n_hw, hw_min, hw_max,  laser_phase ,
                    N_eF, eF_min, eF_max, Tkelvin,eta_smearing,
                    plot_bands, debug_mode,
                    do_gauge_trafo, R_vect_float,
                    do_write_velo, do_write_mep_bands,
                    do_mep, do_kubo, do_ahc, do_opt, do_gyro,
                    verbose=False
                ):

    root_dir, w90_dir   =   dir_setup(  base_dir,   sub_dir)

    #
    print_v(" ^^^^^^^^^^ FeMn TB model setup ^^^^^^^^^^^^^^^^^^",verbose)
    print_v(" *\n *\n *\n",verbose)
    #
    #-------------------------------------------------#
    # TB PARAMETERS                                   #
    #-------------------------------------------------#
    t          = t           # nearest hopping
    strain     = strain     # strain = d0'/d0, strain scale along 111 direction
    J          = J_ex         # local exchange field
    Jz         = Jz       # Zeeman field
    tso        = tso        # spin-orbit coupling
    nns = 1           # only consider the 1st hopping, but it can extend to 2nd, 3rd ...
    #
    d =  strain**2-1    # d>0 (d<0) indicates tensile (compress) strain along 111 direction
    t_intra =  t       # hopping within (111) plane
    #t_inter = t-d;     # hopping between (111) planes
    t_inter =  t/float( strain)**2
    #
    if np.abs( strain-1.)> 1e-3:  # For strained case, the original distance of nearest neighbor ...
        nns = 2                  #   nns = 2        # will change to two different lengths.
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    #
    #
    #-------------------------------------------------#
    # LATTICE                                         #
    #-------------------------------------------------#
    ndim = 3         # 3D system

    #
    print_v('initial lattice :', verbose)
    print_v(latt0, verbose)
    print_v('vol = '+str(vol), verbose)
    print_v("now apply strain="+str(strain)+' ...\n', verbose)
    #
    latt_optimizer               =   FeMn_latt_OPT(strain)
    opt_x, opt_vol,  latt    =   latt_optimizer.optimize_lattice(1e-7)  # solve the strained lattice, by keeping constant volume
    print_v("\n\n\n", verbose)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    #
    #-------------------------------------------------#
    # atomic positions                                #
    #-------------------------------------------------#
    natoms = 4
    num_orb= 2* natoms
    atoms = np.array([  [0.0, 0.0, 0.0],       # 3D system
                        [0.5, 0.5, 0.0],       # four Fe/Mn atoms
                        [0.5, 0.0, 0.5],
                        [0.0, 0.5, 0.5]
                ])
    #
    atoms_cart = np.zeros((natoms,ndim))
    centroid = 0.25*(atoms[0,:]+atoms[1,:]+atoms[2,:]+atoms[3,:]) # the centeroid of a tetrahedron
    centroid_cart = np.zeros(ndim)
    #
    for i in range(natoms):
        for j in range(ndim):
            atoms_cart[i,j] = np.dot(atoms[i,:],latt[:,j])
            centroid_cart[j] = np.dot(centroid[:],latt[:,j])
    #
    print_v('atoms cart:', verbose)
    print_v(atoms_cart, verbose)
    print_v("\n\n\n", verbose)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #
    #-------------------------------------------------#
    # spin configuration                              #
    #-------------------------------------------------#
    if spin_order == '3Q':
        theta   = theta_3q
        phi     = phi_3q
    elif spin_order == '2Q':
        theta   = theta_2q
        phi     = phi_2q
    elif spin_order == '1Q':
        theta   = theta_1q
        phi     = phi_1q
    elif spin_order == '0Q':
        theta   = theta_0q
        phi     = phi_0q
    else:
        theta   = theta_3q
        phi     = phi_3q
        print_v("[write_3q_HR]:   WARNING unknown spin_order ID "+spin_order+" . Will use default (3Q state)", verbose)
    #
    mag = np.zeros((natoms,ndim))
    for i in range(natoms):
        mag[i,0] = np.sin(theta[i])*np.cos(phi[i])
        mag[i,1] = np.sin(theta[i])*np.sin(phi[i])
        mag[i,2] = np.cos(theta[i])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    #
    #
    #-------------------------------------------------#
    # find all atoms in a ncells^ndim supercell       #
    #-------------------------------------------------#
    ncells = 1
    #
    nrpts = (2*ncells+1)**ndim
    totatom = nrpts*natoms
    atomlist = np.zeros((totatom,ndim))
    R = np.zeros((1,ndim))
    index = 0
    #
    # iterate real space neighbours
    for x in range(-ncells, ncells+1,1):     #= -ncells:ncells
        for y in range(-ncells, ncells+1,1):     #= -ncells:ncells
            for z in range(-ncells, ncells+1,1):     #= -ncells:ncells
    #            print(x,' ',y,' ',z)
                for i in range(natoms):
                    R[:] = x*latt[0,:]+y*latt[1,:]+z*latt[2,:]
                    atomlist[index,:] = R + atoms_cart[i,:]
                    index = index + 1
    #
    # find 1st, 2nd, 3rd ... distances.
    index = 0
    aug = np.zeros(natoms*totatom)
    for i in range(natoms):  #= 1:natoms
        for j in range(totatom):  #= 1:natoms
            aug[index] = np.linalg.norm(atoms_cart[i,:]-atomlist[j,:])
            index = index + 1;
    #
    # remove degeneracy
    uni_aug =   np.unique(aug)
    dis = []
    for elem in uni_aug:
        if np.abs(elem)>1e-3:
            dis.append(elem)
    #
    # black magic fuckery
    for i in range(len(dis)):#= 1: length(dis)-1
        if dis[i] > 0:
            for j in range(i+1, len(dis)):
                if  (       np.abs(dis[j]-dis[i])   <   1e-5):
                    dis[j] = 0
    #
    # remove on-site
    dis_nonzero =   []
    for elem in dis:
        if np.abs(elem)>1e-3:
            dis_nonzero.append(elem)
    dis = dis_nonzero
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    #construct H!








    #
    #-------------------------------------------------#
    # WRITE LATTICE                                   #
    #-------------------------------------------------#
    with open(w90_dir+'FeMn_d'+str(strain)+'.latt','w') as out_latt:
        out_latt.write( '#Created on ' +      datetime.date.today().strftime("%d%B%Y")
                        +' [TB FeMn python: '+str(spin_order)+' t='+str(t)+'(eV); strain='+str(strain)
                        +'] \n')
        for i in range(3):
            for j in range(3):
                out_latt.write( '{:.8f}'.format(latt[i,j])+' ')
            out_latt.write('\n')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    #
    #-------------------------------------------------#
    # EXCHANGE TERM H1                                #
    #-------------------------------------------------#
    H1 = np.zeros((num_orb,num_orb))
    for i in range(natoms):                             #= 1:natoms
        pos = np.zeros((natoms,natoms))
        pos[i,i] = 1
        magsig =  mag[i,0]*sigma_x +    \
                  mag[i,1]*sigma_y +    \
                  mag[i,2]*sigma_z
        H1 = H1 + np.kron(magsig,pos)
    H1 = -J*H1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



    #
    #-------------------------------------------------#
    # WRITE WANNIER OUT FILE                          #
    #-------------------------------------------------#
    delta_cart  =   np.zeros(ndim)
    R           =   np.zeros(3)
    R_cart      =   np.zeros(ndim)
    idx_r       =   0


    out_fpath   =   w90_dir + seed_name   + '.dat'
    with open(out_fpath,'w') as out_file:
        #   ^^^^^^^^^^^^^^^^
        #      HEADER
        print_v("[FeMn_WX_setup]: start writing "+out_fpath+" file",verbose)
        out_file.write( 'Created with python on ' +      datetime.date.today().strftime("%d%B%Y")
                        +' [TB FeMn: '+str(spin_order)+' t='+str(t)+'(eV); strain='+str(strain)
                        +'] \n')
        out_file.write('{:11d}\n'.format(num_orb))
        out_file.write('{:11d}\n'.format((2*ncells+1)**ndim))
        k = int(    np.fix(nrpts/15)    )
        for i in range(k):
            for j in range(15):
                out_file.write('{:5d}'.format(1))
            out_file.write('\n')
        #
        for i in range(nrpts-k*15):
            out_file.write('{:5d}'.format(1))
        out_file.write('\n')
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        #
        #
        #   ^^^^^^^^^^^^^^^^
        #       R-SPACE HAMILTONIAN
        print_v("[FeMn_WX_setup]: start writing to body of file "+out_fpath,verbose)
        print_v('\t write R_cell...',verbose)
        for x in range(-ncells, ncells+1,1):    #= -ncells:ncells
            for y in range(-ncells, ncells+1,1):
                for z in range(-ncells, ncells+1,1):
                    R[0] = int(x)
                    R[1] = int(y)
                    R[2] = int(z)
                    idx_r = idx_r +1
                    R_cart[:] = x*latt[0,:] + y*latt[1,:] + z*latt[2,:]
                    print_v("\t\t ... #"+str(idx_r)+': '+str(R_cart),verbose)
                    for i in range(0,num_orb):  #= 1:num_orb
                        for j in range(0,num_orb):
                            HH = 0
                            #
                            # construct H0
                            ii=i
                            jj=j
                            if i>natoms-1:
                                ii = ii-natoms
                            if j>natoms-1:
                                jj = jj-natoms
                            #
                            #print(i,'x',j,' -> ',ii,'x',jj)

                            #-------------------------------------------------#
                            # HOPPING HH                                      #
                            #-------------------------------------------------#
                            if (not(ii ==jj)) and ((i<=(natoms-1) and j<=(natoms-1)) or (i>(natoms-1) and j>(natoms-1))):
                                delta_cart[:] = R_cart[:] + atoms_cart[jj][:] - atoms_cart[ii][:]
                                # first/second nearest-neighbors for without (with) strain
                                if np.linalg.norm(delta_cart)<dis[nns+1]:
                                    #print('delta_cart=',delta_cart)
                                    if np.abs(delta_cart[2])<1e-5:
                                        HH = t_intra
                                    else:
                                        HH = t_inter
                            #
                            # exchange term
                            if x==0 and y==0 and z==0 :
                               HH = HH + H1[i,j]

                #             # SOC term
                #             Hso = zeros(num_orb,num_orb);
                #             if ii~=jj
                #                 delta_cart = R_cart(:)' + atoms_cart(jj,:) - atoms_cart(ii,:);
                #                 if norm(delta_cart)<dis(nns+1)
                #                     pos = zeros(natoms,natoms);
                #                     pos(ii,jj) = 1;
                #                     nsig = n(n_index(ii,jj),1).*sigma_x + ...
                #                            n(n_index(ii,jj),2).*sigma_y + ...
                #                            n(n_index(ii,jj),3).*sigma_z;
                #                     Hso = 1i*tso*fac*v(ii,jj)*kron(nsig,pos);
                #                     HH = HH + Hso(i,j);
                #                 end
                #             end
                           # fprintf(f,'#5i#5i#5i#5i#5i#12.6f#12.6f\n', R(:),j,i,real(HH),imag(HH));
                            line    = '\t' +  str(int(R[0]))+   '\t' +    str(int(R[1]))  +'\t'+       str(int(R[2]))        +'\t'
                            line    =   line    +   str(j+1)+'\t'+str(i+1)                            +'\t'
                            line    =   line    +   '\t{:.6f}'.format(np.real(HH))
                            line    =   line    +   '\t{:.6f}'.format(np.imag(HH))
                            #
                            line    =   '{:5d}{:5d}{:5d}{:5d}{:5d}{:12.6f}{:12.6f}\n'.format(int(R[0]),int(R[1]),int(R[2]),j+1,i+1,np.real(HH),np.imag(HH))
                            #
                            out_file.write(line)
        out_file.close()
        print_v('wrote '+str(idx_r)+" cells to Hamiltonian",verbose)

        #
        #-------------------------------------------------#
        # WRITE WANNIER INTERPOLATION INPUT               #
        #-------------------------------------------------#


        ax  =   latt[0,:]
        ay  =   latt[1,:]
        az  =   latt[2,:]
        a0  =   1

        file_path   =   root_dir
        #def write_mepInterp_input(  file_path,valence_bands, ax, ay, az, a0, mp_grid, seed_name,
        #                    kubo_tol, n_hw, hw_min, hw_max,  laser_phase ,N_eF, eF_min, eF_max, Tkelvin,eta_smearing,
        #                    plot_bands, debug_mode, do_gauge_trafo, R_vect_float    , do_write_velo,    do_write_mep_bands,
        #                    do_mep, do_kubo, do_ahc, do_opt, do_gyro
        #                ):

        write_mepInterp_input(
                            file_path,valence_bands, ax, ay, az, a0, mp_grid, seed_name,
                            kubo_tol, n_hw, hw_min, hw_max,  laser_phase ,N_eF, eF_min, eF_max, Tkelvin,eta_smearing,
                            plot_bands, debug_mode, do_gauge_trafo, R_vect_float    , do_write_velo,    do_write_mep_bands,
                            do_mep, do_kubo, do_ahc, do_opt, do_gyro, verbose
                        )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#



def loop_strain(dmin, dmax, n_d):
    base_dir    =   'newRun_'+datetime.date.today().strftime("%d%B%Y")
    # 
    print('setup FeMn model at different strain values')   
    for strain in np.linspace(dmin,dmax,n_d):
        #
        print("\traw strain:\t",strain)
        #   truncate everything after 3 digit
        strain  =   strain * 1000
        strain  =   np.trunc(strain)
        strain  =   strain / 1000.
        #
        #   name of subdir where calc. will be perfomed in
        sub_dir =   'd{:.3f}'.format(strain)
        #
        write_3q_HR(        base_dir    =   base_dir           ,
                            sub_dir     =   sub_dir             ,
                            seed_name   =    'wf1_hr'           ,
                            t           =   -    1.0            ,
                            strain      =       strain           ,
                            J_ex        =   -    1.0             ,
                            Jz          =         0             ,
                            tso         =         0             ,
                            spin_order  =       '3Q'            ,
                            #
                            mp_grid     =  [200,200,200]         ,
                            kubo_tol    =           1e-5         ,
                            valence_bands=          2            ,
                            #
                            n_hw=                   120           ,
                            hw_min=                 0             ,
                            hw_max=                 6             ,
                            laser_phase=            1             ,
                            #
                            N_eF=                   1             ,
                            eF_min=                 0             ,
                            eF_max=                 0             ,
                            Tkelvin=                0             ,
                            eta_smearing=           0.1           ,
                            #
                            plot_bands=              False         ,
                            debug_mode=              False         ,
                            do_gauge_trafo=          True          ,
                            R_vect_float=            False         ,
                            do_write_velo=           False         ,
                            do_write_mep_bands=      True          ,
                            #
                            do_mep=             True                ,
                            do_kubo=            False                ,
                            do_ahc=             True                ,
                            do_opt=             False               ,
                            do_gyro=            False               ,
                            verbose=            False
                    )
        #
        print("\tfinished setup of folder:\t\t->",sub_dir)
    print('\nfinished setting up ./'+base_dir)




loop_strain(0.9,1.1,5)












