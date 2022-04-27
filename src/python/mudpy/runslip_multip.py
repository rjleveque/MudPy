'''
Diego Melgar, 01/2014
Runtime file for forward modeling and inverse kinematic slip inversions
'''

# Need environment to pass into Popen with paths set for 
# importing mudpy and running fk.pl
import os,sys
MUD = os.environ['MUD']
fkdir = '%s/src/fk/' % MUD
mud_env = os.environ.copy()
mud_env['PATH'] = mud_env['PATH'] + ':' + fkdir
mud_env['PYTHONPATH'] = mud_env['PYTHONPATH'] + ':' + '%s/src/python' % MUD
print('+++ Will call Popen with PATH = ',mud_env['PATH'])
print('+++ and PYTHONPATH = ',mud_env['PYTHONPATH'])


#Initalize project folders
def init(home,project_name):
    '''
    Initalizes file structure for a new problem
    
    IN:
        home: What dir will you be working from
        project_name: What name will you give this problem
        
    OUT:
        Nothing
    '''
    from shutil import rmtree,copy
    from os import path,makedirs,environ
    clob='y'
    proj_dir=home+project_name+'/'
    if path.exists(proj_dir):  #Path exists, clobber?
        clob=input('Project directory exists, clobber (y/n)?')
        if clob =='y' or clob == 'Y': #Clobber baby
            clob=input('This will delete everything in this project directory, so, take a minute, think about it: clobber (y/n)?')
            if clob == 'y' or clob == 'Y':
                rmtree(proj_dir)
            else: #Leave direcory alone
                print('Phew, almost shot yourself in the foot there didn\'t you?')
        else: #Leave directory alone
            print('Phew, almost shot yourself in the foot there didn\'t you?')
    if clob == 'y' or clob == 'Y':
        makedirs(proj_dir)
        #And make the subdirectories
        makedirs(proj_dir+'GFs')
        makedirs(proj_dir+'GFs/static')
        makedirs(proj_dir+'GFs/dynamic')
        makedirs(proj_dir+'GFs/matrices')
        makedirs(proj_dir+'GFs/tsunami')
        makedirs(proj_dir+'GFs/STFs')
        makedirs(proj_dir+'data/waveforms')
        makedirs(proj_dir+'data/statics')
        makedirs(proj_dir+'data/station_info')
        makedirs(proj_dir+'data/model_info')
        makedirs(proj_dir+'structure')
        makedirs(proj_dir+'plots')
        makedirs(proj_dir+'scripts')
        makedirs(proj_dir+'forward_models')
        makedirs(proj_dir+'output/inverse_models')
        makedirs(proj_dir+'output/inverse_models/statics')
        makedirs(proj_dir+'output/inverse_models/waveforms')
        makedirs(proj_dir+'output/inverse_models/models')
        makedirs(proj_dir+'output/forward_models')
        makedirs(proj_dir+'logs')
        makedirs(proj_dir+'analysis')
        makedirs(proj_dir+'analysis/frequency')
        #Copy templates into appropriate files
        mudpy=environ['MUD']+'/run/'
        copy(mudpy+'template.fault',proj_dir+'data/model_info/')
        copy(mudpy+'template.gflist',proj_dir+'data/station_info/')
        copy(mudpy+'template.sta',proj_dir+'data/station_info/')
        copy(mudpy+'template.mod',proj_dir+'structure/')


#Extract fault geometry from rupture file
def rupt2fault(home,project_name,rupture_name):
    '''
    Make fault file from user provided forward model rupture file
    '''
    from numpy import loadtxt,savetxt,c_
    
    print('Assembling fault file from rupture file')
    rupt=loadtxt(home+project_name+'/forward_models/'+rupture_name,ndmin=2)
    fault=c_[rupt[:,0],rupt[:,1],rupt[:,2],rupt[:,3],rupt[:,4],rupt[:,5],rupt[:,6],rupt[:,7],rupt[:,10],rupt[:,11]]
    savetxt(home+project_name+'/data/model_info/'+rupture_name.split('.')[0]+'.fault', \
            fault,fmt='%i\t%.6f\t%.6f\t%.6f\t%.2f\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f')





# Run green functions          
def make_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static,tsunami,
            hot_start,dk,pmin,pmax,kmax,okada=False):
    '''
    [Original serial version - deprecated.]
    This routine set's up the computation of GFs for each subfault to all stations.
    The GFs are impulse sources, they don't yet depend on strike and dip.
    
    IN:
        home: Home directory
        project_name: Name fo the problem
        station_file: File with coordinates of stations
        fault_name: Name of fault file
        model_Name: Name of Earth structure model file
        dt: Desired sampling itnerval for waveforms
        NFFT: No. of samples requested in waveform (must be power of 2)
        static: =0 if computing full waveforms, =1 if computing only the static field
        hot_start: =k if you want to start computations at k-th subfault, =0 to compute all

        
    OUT:
        Nothing
    '''
    import time,glob
    from mudpy import green
    from numpy import loadtxt
    from shutil import rmtree,copy
    from os import chdir,path,makedirs,remove
    import datetime
    import gc
    
    if okada==False:
        tic=time.time()
        model_path=home+project_name+'/structure/'
        green_path=home+project_name+'/GFs/'
        station_file=home+project_name+'/data/station_info/'+station_file 
        fault_file=home+project_name+'/data/model_info/'+fault_name  
        logpath=home+project_name+'/logs/'
        #log time
        now=datetime.datetime.now()
        now=now.strftime('%b-%d-%H%M')
        chdir(model_path)
        #Load source model for station-event distance computations
        source=loadtxt(fault_file,ndmin=2)
        for k in range(hot_start,source.shape[0]):
            #Run comptuation for 1 subfault
            log=green.run_green(source[k,:],station_file,model_name,dt,NFFT,static,dk,pmin,pmax,kmax)
            #Write log
            f=open(logpath+'make_green.'+now+'.log','a')    
            f.write(log)
            f.close()
            #Move to correct directory
            strdepth='%.4f' % source[k,3]
            subfault=str(int(source[k,0])).rjust(4,'0') #subfault=rjust(str(int(source[k,0])),4,'0')
            if static==0 and tsunami==False:
                #Move results to dynamic GF dir
                dirs=glob.glob('*.mod_'+strdepth)
                #Where am I writting this junk too?
                outgreen=green_path+'/dynamic/'+path.split(dirs[0])[1]+'.sub'+subfault
                #Check if GF subdir already exists
                if path.exists(outgreen)==False:
                    #It doesn't, make it, don't be lazy
                    makedirs(outgreen)
                #Now copy GFs in, this will OVERWRITE EXISTING GFs of the same name
                flist=glob.glob(dirs[0]+'/*')
                for k in range(len(flist)):
                    copy(flist[k],outgreen)
                #Cleanup
                rmtree(dirs[0])
                gc.collect()
            elif static==0 and tsunami==True: #Tsunami GFs
                #Move results to tsunami GF dir
                dirs=glob.glob('*.mod_'+strdepth)
                #Where am I writting this junk too?
                outgreen=green_path+'/tsunami/'+path.split(dirs[0])[1]+'.sub'+subfault
                #Check if GF subdir already exists
                if path.exists(outgreen)==False:
                    #It doesn't, make it, don't be lazy
                    makedirs(outgreen)
                #Now copy GFs in, this will OVERWRITE EXISTING GFs of the same name
                flist=glob.glob(dirs[0]+'/*')
                for k in range(len(flist)):
                    copy(flist[k],outgreen)
                #Cleanup
                rmtree(dirs[0])
                gc.collect()
            else:  #Static GFs
                copy('staticgf',green_path+'static/'+model_name+'.static.'+strdepth+'.sub'+subfault)
                #Cleanup
                remove('staticgf')     
        #How long was I working for?
        toc=time.time()
        print('GFs computed in '+str((toc-tic)/60)+' minutes...')
    else:
        print('GFs not necessary when using an elastic halfspace, exiting make_green')


def make_green_one_subfault(ksource,home,project_name,station_file,
                            model_name,dt,NFFT,static,dk,pmin,pmax,
                            kmax,tsunami,insar,source):
                            
    """
    Make the Greens function for one subfault defined by source[ksource,:].
    """                 
    import subprocess
    from os import chdir
    from shutil import copy,rmtree
    from numpy import genfromtxt,zeros
    from shlex import split
    from shutil import copy
    from glob import glob
    from mudpy.green import src2sta
    import os
    from multiprocessing import current_process
    
    #print('+++ in make_green_one_subfault, ksource = %i' % ksource)
    ppool = current_process()
    print('Pool  Process %i now working on subfault %i' \
            % (ppool.pid, source[ksource,0]))
    
    #Where should I be working boss?
    depth='%.4f' % source[ksource,3]
    subfault=str(int(source[ksource,0])).rjust(4,'0')
    if tsunami==False and static==0:
        subfault_folder=home+project_name+'/GFs/dynamic/'+model_name+'_'+depth+'.sub'+subfault
    elif tsunami==True and static==1:
        subfault_folder=home+project_name+'/GFs/tsunami/'+model_name+'_'+depth+'.sub'+subfault
    elif static==1:
        subfault_folder=home+project_name+'/GFs/static/'+model_name+'_'+depth+'.sub'+subfault
    
    #Check if subfault folder exists, if not create it
    if os.path.exists(subfault_folder+'/')==False:
        os.makedirs(subfault_folder+'/')
    
    #Copy velocity model file
    copy(home+project_name+'/structure/'+model_name,subfault_folder+'/'+model_name)
    #Move to work folder
    chdir(subfault_folder)
    #Get station distances to source
    d,az=src2sta(station_file,source[ksource,:])
    #Make distance string for system call
    diststr=''
    for k in range(len(d)):
        diststr=diststr+' %.6f' % d[k] #Truncate distance to 6 decimal palces (meters)
    # Keep the user informed, lest he get nervous
    #print('MPI: processor #',rank,'is now working on subfault',int(source[ksource,0]),'(',ksource+1,'/',len(source),')')
    

    #Make the calculation
    if static==0: #Compute full waveform
        command=split("fk.pl -M"+model_name+"/"+depth+"/f -N"+str(NFFT) \
                    +"/"+str(dt)+'/1/'+repr(dk)+' -P'+repr(pmin)+'/' \
                    +repr(pmax)+'/'+repr(kmax)+diststr)
        p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                           env=mud_env)
        #print('Popen Process %i will run fk on subfault %i' \
        #        % (p.pid,source[ksource,0]))
        p.communicate() 
        # Move files up one level and delete folder created by fk
        files_list=glob(subfault_folder+'/'+model_name+'_'+depth+'/*.grn*')
        for f in files_list:
            newf=subfault_folder+'/'+f.split('/')[-1]
            #print('+++ copying %s %s' % (f,newf))
            copy(f,newf)
        rmtree(subfault_folder+'/'+model_name+'_'+depth)
    else: #Compute only statics
        if insar==True:
            suffix='insar'
        else:
            suffix='gps'
        write_file=subfault_folder+'/'+model_name+'.static.'+depth+'.sub'+subfault+'.'+suffix
        command=split("fk.pl -M"+model_name+"/"+depth+"/f -N1 "+diststr)
        file_is_empty=True
        while file_is_empty:
            p=subprocess.Popen(command,stdout=open(write_file,'w'),stderr=subprocess.PIPE)
            p.communicate()
            if os.stat(write_file).st_size!=0: #File is NOT empty
                file_is_empty=False
            else:
                print('Warning: I just had a mini-seizure and made an empty GF file on first try, re-running')
        #If file is empty run again   
         
    #print('Popen Process %i done running fk on subfault %i' \
    #        % (p.pid,source[ksource,0]))
            
    ppool = current_process()
    print('Pool  Process %i done working on subfault %i' \
            % (ppool.pid, source[ksource,0]))
    # end of make_green_one_subfault



def make_parallel_green_multip(home,project_name,station_file,fault_name,
            model_name,dt,NFFT,static,tsunami,
            hot_start,dk,pmin,pmax,kmax,ncpus,insar=False,okada=False):
    '''
    [New parallel version using multiprocessing.pool to farm out subfaults.]
    
    This routine set's up the computation of GFs for each subfault to all stations.
    The GFs are impulse sources, they don't yet depend on strike and dip.
    
    IN:
        home: Home directory
        project_name: Name fo the problem
        station_file: File with coordinates of stations
        fault_name: Name of fault file
        model_Name: Name of Earth structure model file
        dt: Desired sampling itnerval for waveforms
        NFFT: No. of samples requested in waveform (must be power of 2)
        static: =0 if computing full waveforms, =1 if computing only the static field
        hot_start: =k if you want to start computations at k-th subfault, =0 to compute all

        
    OUT:
        Nothing
    '''
    
    from numpy import loadtxt,arange,savetxt
    from os import path,makedirs,environ
    from shlex import split
    import subprocess
    from multiprocessing import Pool, current_process
    
    green_path=home+project_name+'/GFs/'
    station_file=home+project_name+'/data/station_info/'+station_file 
    fault_file=home+project_name+'/data/model_info/'+fault_name  
    #Load source model for station-event distance computations
    source=loadtxt(fault_file,ndmin=2)
    
    #Create all output folders
    for k in range(hot_start,source.shape[0]):
        strdepth='%.4f' % source[k,3]
        subfault=str(k+1).rjust(4,'0')
        if static==0 and tsunami==False:
            subfault_folder=green_path+'dynamic/'+model_name+'_'+strdepth+'.sub'+subfault
            if path.exists(subfault_folder)==False:
                #It doesn't, make it, don't be lazy
                makedirs(subfault_folder)               
#        elif static==0 and tsunami==True: #Tsunami GFs
#            subfault_folder=green_path+'tsunami/'+model_name+'_'+strdepth+'.sub'+subfault
#            if path.exists(subfault_folder)==False:
#                #It doesn't, make it, don't be lazy
#                makedirs(subfault_folder)
        elif static==1 and tsunami==True: #Tsunami GFs
            subfault_folder=green_path+'tsunami/'+model_name+'_'+strdepth+'.sub'+subfault
            if path.exists(subfault_folder)==False:
                #It doesn't, make it, don't be lazy
                makedirs(subfault_folder)


    # Adapted from parallel.run_parallel_green:
    #What parameters are we using?
    if 1:
        out='''Running all processes with:
        home = %s
        project_name = %s
        station_file = %s
        model_name = %s
        static = %s
        tsunami = %s
        dt = %.3f
        NFFT = %d
        dk = %.3f
        pmin = %.3f
        pmax = %.3f
        kmax = %.3f
        insar = %s
        ''' %(home,project_name,station_file,model_name,str(static),\
              str(tsunami),dt,NFFT,dk,pmin,pmax,kmax,str(insar))
        print(out)

    # Make a list of arguments needed for make_green_one_subfault
    # Only the first argument ksource changes as we loop through all the 
    # subfaults. (Which will be done in parallel with multiprocessing.Pool.)
    
    ksources = list(range(hot_start,source.shape[0]))
    print('Will make Greens functions for %i subfaults...' % len(ksources))
    
    all_args = []
    for ksource in ksources:
        args = (ksource,home,project_name,station_file,
                               model_name,dt,NFFT,static,dk,pmin,pmax,
                               kmax,tsunami,insar,source)
        all_args.append(args)
        
    # now all_args is a list and each element of the list is a tuple with all
    # the arguments needed for make_green_one_subfault.  
    # Farm these out to the requested number of threads, ncpus:
    
    with Pool(processes=ncpus) as pool:
        pool.starmap(make_green_one_subfault, all_args)
    
        
        
        
def make_parallel_teleseismics_green(home,project_name,station_file,fault_name,model_name,teleseismic_vel_mod,time_epi,endtime,ncpus,hot_start=0):
    '''
    [mpi version, not yet rewritten for multip.]
    This routine set's up the computation of GFs for each subfault to all stations.
    The GFs are impulse sources, they don't yet depend on strike and dip.
    
    IN:
        home: Home directory
        project_name: Name fo the problem
        station_file: File with coordinates of stations
        fault_name: Name of fault file
        model_Name: Name of Earth structure model file
        dt: Desired sampling itnerval for waveforms
        NFFT: No. of samples requested in waveform (must be power of 2)
        static: =0 if computing full waveforms, =1 if computing only the static field
        hot_start: =k if you want to start computations at k-th subfault, =0 to compute all

        
    OUT:
        Nothing
    '''
    from numpy import arange,savetxt,genfromtxt
    from os import path,makedirs,environ
    from shlex import split
    import subprocess
    
    
    green_path=home+project_name+'/GFs/'
    station_file=home+project_name+'/data/station_info/'+station_file 
    fault_file=home+project_name+'/data/model_info/'+fault_name  
    #Load source model for station-event distance computations
    source=genfromtxt(fault_file)
   
    #Create all output folders
    for k in range(len(source)):
        
        strdepth='%.4f' % source[k,3]
        subfault=str(k+1).rjust(4,'0')
        
        subfault_folder=green_path+'dynamic/'+model_name+'_'+strdepth+'.sub'+subfault
        if path.exists(subfault_folder)==False:
            #It doesn't, make it, don't be lazy
            makedirs(subfault_folder)               

    #Create individual source files
    for k in range(ncpus):
        i=arange(k+hot_start,len(source),ncpus)
        mpi_source=source[i,:]
        fmt='%d\t%10.6f\t%10.6f\t%.8f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f'
        savetxt(home+project_name+'/data/model_info/mpi_source.'+str(k)+'.fault',mpi_source,fmt=fmt)
    
    #Make mpi system call
    print("MPI: Starting GFs computation on", ncpus, "CPUs\n")
    mud_source=environ['MUD']+'/src/python/mudpy/'


    mpi='mpiexec -n '+str(ncpus)+' python '+mud_source+'parallel.py run_parallel_teleseismics_green '+home+' '+project_name+' '+str(time_epi)+' '+station_file+' '+model_name+' '+teleseismic_vel_mod+' '+str(endtime)
    print(mpi)
    mpi=split(mpi)
    p=subprocess.Popen(mpi, env=mud_env)
    p.communicate()
        
        


   

#Now make synthetics for source/station pairs
def make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,tsunami,beta,
                    hot_start,time_epi,impulse=False,okada=False,mu=45e9,insar=False):
    '''
    [Original serial version -- deprecated.]
    This routine will take the impulse response (GFs) and pass it into the routine that will
    convovle them with the source time function according to each subfaults strike and dip.
    The result fo this computation is a time series dubbed a "synthetic"
    
    IN:
        home: Home directory
        project_name: Name fo the problem
        station_file: File with coordinates of stations
        fault_name: Name of fault file
        model_Name: Name of Earth structure model file
        integrate: =0 if you want output to be velocity, =1 if you want output to be displacement
        static: =0 if computing full waveforms, =1 if computing only the static field
        hot_start: =k if you want to start computations at k-th subfault, =0 to compute all

        
    OUT:
        Nothing
    '''
    from mudpy import green
    import datetime
    from numpy import loadtxt
    import gc
    
    green_path=home+project_name+'/GFs/'
    station_file=home+project_name+'/data/station_info/'+station_file
    fault_file=home+project_name+'/data/model_info/'+fault_name
    logpath=home+project_name+'/logs/'
    #Time for log file
    now=datetime.datetime.now()
    now=now.strftime('%b-%d-%H%M')
    #First read fault model file
    source=loadtxt(fault_file,ndmin=2)
    #Now compute synthetics please, one sub fault at a time
    for k in range(hot_start,source.shape[0]):
        print('ksource = ' + str(k))
        subfault=str(k+1).rjust(4,'0')
        log=green.run_syn(home,project_name,source[k,:],station_file,green_path,model_name,integrate,static,tsunami,
                subfault,time_epi,beta,impulse,okada,mu,insar=insar)
        f=open(logpath+'make_synth.'+now+'.log','a')
        f.write(log)
        f.close()
        gc.collect()
        
        
#Now make synthetics for source/station pairs
def make_parallel_synthetics(home,project_name,station_file,
            fault_name,model_name,integrate,static,tsunami,beta,
            hot_start,time_epi,ncpus,custom_stf,impulse=False,
            insar=False,okada=False,mu=45e9):
    '''
    [Old mpi version]
    This routine will take the impulse response (GFs) and pass it into the routine that will
    convovle them with the source time function according to each subfaults strike and dip.
    The result fo this computation is a time series dubbed a "synthetic"
    
    IN:
        home: Home directory
        project_name: Name fo the problem
        station_file: File with coordinates of stations
        fault_name: Name of fault file
        model_Name: Name of Earth structure model file
        integrate: =0 if you want output to be velocity, =1 if you want output to be displacement
        static: =0 if computing full waveforms, =1 if computing only the static field
        hot_start: =k if you want to start computations at k-th subfault, =0 to compute all
        coord_type: =0 if problem is in cartesian coordinates, =1 if problem is in lat/lon
        
    OUT:
        Nothing
    '''
    from numpy import arange,savetxt
    import datetime
    from numpy import loadtxt
    import subprocess
    from shlex import split
    from os import environ
    
    station_file=home+project_name+'/data/station_info/'+station_file
    fault_file=home+project_name+'/data/model_info/'+fault_name
    #Time for log file
    now=datetime.datetime.now()
    now=now.strftime('%b-%d-%H%M')
    #First read fault model file
    source=loadtxt(fault_file,ndmin=2)
    #Create individual source files
    for k in range(ncpus):
        i=arange(k+hot_start,len(source),ncpus)
        mpi_source=source[i,:]
        fmt='%d\t%10.6f\t%10.6f\t%.8f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f'
        savetxt(home+project_name+'/data/model_info/mpi_source.'+str(k)+'.fault',mpi_source,fmt=fmt)
    #Make mpi system call
    print("MPI: Starting synthetics computation on", ncpus, "CPUs\n")
    mud_source=environ['MUD']+'/src/python/mudpy/'
    mpi='mpiexec -n '+str(ncpus)+' python '+mud_source+'parallel.py run_parallel_synthetics '+home+' '+project_name+' '+station_file+' '+model_name+' '+str(integrate)+' '+str(static)+' '+str(tsunami)+' '+str(time_epi)+' '+str(beta)+' '+str(custom_stf)+' '+str(impulse)+' '+str(insar)+' '+str(okada)+' '+str(mu)
    print(mpi)
    mpi=split(mpi)
    p=subprocess.Popen(mpi, env=mud_env)
    p.communicate()
        
        
def make_synthetics_one_subfault(one_source,home,project_name,station_file,
                    fault_name,model_name,integrate,static,tsunami,beta,
                    hot_start,time_epi,ncpus,custom_stf,impulse=False,
                    insar=False,okada=False,mu_okada=45e9):
                    
    import os
    import subprocess
    from pandas import DataFrame as df
    from mudpy.forward import get_mu
    from numpy import array,genfromtxt,loadtxt,savetxt,log10,zeros,sin,cos,ones,deg2rad
    from obspy import read
    from shlex import split
    from mudpy.green import src2sta,rt2ne,origin_time,okada_synthetics
    from glob import glob
    from mudpy.green import silentremove
    from os import remove
    from multiprocessing import current_process
    

    #Constant parameters
    rakeDS=90+beta #90 is thrust, -90 is normal
    rakeSS=0+beta #0 is left lateral, 180 is right lateral
    tb=50 #Number of samples before first arrival (should be 50, NEVER CHANGE, if you do then adjust in fk.pl)
    
    #Figure out custom STF
    if custom_stf:
        if custom_stf.lower()!='none':
            custom_stf=home+project_name+'/GFs/STFs/'+custom_stf
        else:
            custom_stf=None
    
    #Load structure
    model_file=home+project_name+'/structure/'+model_name
    structure=loadtxt(model_file,ndmin=2)
    
    #this keeps track of statics dataframe
    write_df=False
        
    #source=mpi_source[ksource,:]
    # in this version, one_source was passed in:
    source = one_source
    
    #Parse the source information
    num=str(int(source[0])).rjust(4,'0')
    xs=source[1]
    ys=source[2]
    zs=source[3]
    strike=source[4]
    dip=source[5]
    rise=source[6]
    if impulse==True:
        duration=0
    else:
        duration=source[7]
    ss_length=source[8]
    ds_length=source[9]
    ss_length_in_km=ss_length/1000.
    ds_length_in_km=ds_length/1000.
    strdepth='%.4f' % zs
    subfault=str(int(source[0])).rjust(4,'0')
    if static==0 and tsunami==0:  #Where to save dynamic waveforms
        green_path=home+project_name+'/GFs/dynamic/'+model_name+"_"+strdepth+".sub"+subfault+"/"
    if static==1 and tsunami==1:  #Where to save dynamic waveforms
        green_path=home+project_name+'/GFs/tsunami/'+model_name+"_"+strdepth+".sub"+subfault+"/"
    if static==1 and tsunami==0:  #Where to save statics
        green_path=home+project_name+'/GFs/static/'+model_name+"_"+strdepth+".sub"+subfault+"/"
    staname=genfromtxt(station_file,dtype="U",usecols=0)
    if staname.shape==(): #Single staiton file
        staname=array([staname])
    #Compute distances and azimuths
    d,az,lon_sta,lat_sta=src2sta(station_file,source,output_coordinates=True)
    
    #Get moment corresponding to 1 meter of slip on subfault
    mu=get_mu(structure,zs)
    Mo=mu*ss_length*ds_length*1.0
    Mw=(2./3)*(log10(Mo)-9.1)
    
    #Move to output folder
    os.chdir(green_path)
    #print('Processor '+str(rank)+' is working on subfault '+str(int(source[0]))+' and '+str(len(d))+' stations ')
    
    ppool = current_process()
    print('Pool  Process %i creating synthetics for subfault %i at %i stations'\
            % (ppool.pid, one_source[0],len(d)))   
                
                
    #This is looping over "sites"
    for k in range(len(d)):
        
        if static==0: #Compute full waveforms
            diststr='%.6f' % d[k] #Need current distance in string form for external call
            #Form the strings to be used for the system calls according to user desired options
            if integrate==1: #Make displ.
                #First Stike-Slip GFs
                if custom_stf==None:
                    commandSS="syn -I -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeSS)+" -D"+str(duration)+ \
                        "/"+str(rise)+" -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".SS.disp.x -G"+green_path+diststr+".grn.0"
                    commandSS=split(commandSS) #Split string into lexical components for system call
                    #Now dip slip
                    commandDS="syn -I -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeDS)+" -D"+str(duration)+ \
                        "/"+str(rise)+" -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".DS.disp.x -G"+green_path+diststr+".grn.0"
                    commandDS=split(commandDS)
                else:
                    commandSS="syn -I -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeSS)+" -S"+custom_stf+ \
                        " -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".SS.disp.x -G"+green_path+diststr+".grn.0"
                    commandSS=split(commandSS) #Split string into lexical components for system call
                    #Now dip slip
                    commandDS="syn -I -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeDS)+" -S"+custom_stf+ \
                        " -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".DS.disp.x -G"+green_path+diststr+".grn.0"
                    commandDS=split(commandDS)
            else: #Make vel.
                #First Stike-Slip GFs
                if custom_stf==None:
                    commandSS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeSS)+" -D"+str(duration)+ \
                        "/"+str(rise)+" -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".SS.vel.x -G"+green_path+diststr+".grn.0"
                    commandSS=split(commandSS)
                    #Now dip slip
                    commandDS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeDS)+" -D"+str(duration)+ \
                        "/"+str(rise)+" -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".DS.vel.x -G"+green_path+diststr+".grn.0"
                    commandDS=split(commandDS)
                else:
                    commandSS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeSS)+" -S"+custom_stf+ \
                        " -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".SS.vel.x -G"+green_path+diststr+".grn.0"
                    commandSS=split(commandSS)
                    #Now dip slip
                    commandDS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeDS)+" -S"+custom_stf+ \
                        " -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".DS.vel.x -G"+green_path+diststr+".grn.0"
                    commandDS=split(commandDS)
            #Run the strike- and dip-slip commands (make system calls)
            p=subprocess.Popen(commandSS, env=mud_env)
            p.communicate() 
            # print('Popen Process %i will run syn SS on subfault %i' \
            #                 % (p.pid,source[0]))
            p=subprocess.Popen(commandDS, env=mud_env)
            p.communicate()
            # print('Popen Process %i will run syn DS on subfault %i' \
            #                 % (p.pid,source[0]))
            #Result is in RTZ system (+Z is down) 
            #rotate to NEZ with +Z up and scale to m or m/s
            if integrate==1: #'tis displacememnt
                #Strike slip
                if duration>0: #Is there a source time fucntion? Yes!
                    r=read(staname[k]+".subfault"+num+'.SS.disp.r')
                    t=read(staname[k]+".subfault"+num+'.SS.disp.t')
                    z=read(staname[k]+".subfault"+num+'.SS.disp.z')
                else: #No! This is the impulse response!
                    r=read(staname[k]+".subfault"+num+'.SS.disp.ri')
                    t=read(staname[k]+".subfault"+num+'.SS.disp.ti')
                    z=read(staname[k]+".subfault"+num+'.SS.disp.zi')
                ntemp,etemp=rt2ne(r[0].data,t[0].data,az[k])
                #Scale to m and overwrite with rotated waveforms
                n=r.copy()
                n[0].data=ntemp/100
                e=t.copy()
                e[0].data=etemp/100
                z[0].data=z[0].data/100
                # get rid of numerical "noise" in the first tb samples
                n[0].data[0:tb]=0
                e[0].data[0:tb]=0
                z[0].data[0:tb]=0
                n=origin_time(n,time_epi,tb)
                e=origin_time(e,time_epi,tb)
                z=origin_time(z,time_epi,tb)
                n.write(staname[k]+".subfault"+num+'.SS.disp.n',format='SAC')
                e.write(staname[k]+".subfault"+num+'.SS.disp.e',format='SAC')
                z.write(staname[k]+".subfault"+num+'.SS.disp.z',format='SAC')
                silentremove(staname[k]+".subfault"+num+'.SS.disp.r')
                silentremove(staname[k]+".subfault"+num+'.SS.disp.t')
                if impulse==True:
                    silentremove(staname[k]+".subfault"+num+'.SS.disp.ri')
                    silentremove(staname[k]+".subfault"+num+'.SS.disp.ti')
                    silentremove(staname[k]+".subfault"+num+'.SS.disp.zi')
                #Dip Slip
                if duration>0: #Is there a source time fucntion? Yes!
                    r=read(staname[k]+".subfault"+num+'.DS.disp.r')
                    t=read(staname[k]+".subfault"+num+'.DS.disp.t')
                    z=read(staname[k]+".subfault"+num+'.DS.disp.z')
                else: #No! This is the impulse response!
                    r=read(staname[k]+".subfault"+num+'.DS.disp.ri')
                    t=read(staname[k]+".subfault"+num+'.DS.disp.ti')
                    z=read(staname[k]+".subfault"+num+'.DS.disp.zi')
                ntemp,etemp=rt2ne(r[0].data,t[0].data,az[k])
                n=r.copy()
                n[0].data=ntemp/100
                e=t.copy()
                e[0].data=etemp/100
                z[0].data=z[0].data/100
                n=origin_time(n,time_epi,tb)
                e=origin_time(e,time_epi,tb)
                z=origin_time(z,time_epi,tb)
                n.write(staname[k]+".subfault"+num+'.DS.disp.n',format='SAC')
                e.write(staname[k]+".subfault"+num+'.DS.disp.e',format='SAC')
                z.write(staname[k]+".subfault"+num+'.DS.disp.z',format='SAC')
                silentremove(staname[k]+".subfault"+num+'.DS.disp.r')
                silentremove(staname[k]+".subfault"+num+'.DS.disp.t')
                if impulse==True:
                    silentremove(staname[k]+".subfault"+num+'.DS.disp.ri')
                    silentremove(staname[k]+".subfault"+num+'.DS.disp.ti')
                    silentremove(staname[k]+".subfault"+num+'.DS.disp.zi')
            else: #Waveforms are velocity, as before, rotate from RT-Z to NE+Z and scale to m/s
                #Strike slip
                if duration>0: #Is there a source time fucntion? Yes!
                    r=read(staname[k]+".subfault"+num+'.SS.vel.r')
                    t=read(staname[k]+".subfault"+num+'.SS.vel.t')
                    z=read(staname[k]+".subfault"+num+'.SS.vel.z')
                else: #No! This is the impulse response!
                    r=read(staname[k]+".subfault"+num+'.SS.vel.ri')
                    t=read(staname[k]+".subfault"+num+'.SS.vel.ti')
                    z=read(staname[k]+".subfault"+num+'.SS.vel.zi')
                ntemp,etemp=rt2ne(r[0].data,t[0].data,az[k])
                n=r.copy()
                n[0].data=ntemp/100
                e=t.copy()
                e[0].data=etemp/100
                z[0].data=z[0].data/100
                n=origin_time(n,time_epi,tb)
                e=origin_time(e,time_epi,tb)
                z=origin_time(z,time_epi,tb)
                n.write(staname[k]+".subfault"+num+'.SS.vel.n',format='SAC')
                e.write(staname[k]+".subfault"+num+'.SS.vel.e',format='SAC')
                z.write(staname[k]+".subfault"+num+'.SS.vel.z',format='SAC')
                silentremove(staname[k]+".subfault"+num+'.SS.vel.r')
                silentremove(staname[k]+".subfault"+num+'.SS.vel.t')
                if impulse==True:
                    silentremove(staname[k]+".subfault"+num+'.SS.vel.ri')
                    silentremove(staname[k]+".subfault"+num+'.SS.vel.ti')
                    silentremove(staname[k]+".subfault"+num+'.SS.vel.zi')
                #Dip Slip
                if duration>0: #Is there a source time fucntion? Yes!
                    r=read(staname[k]+".subfault"+num+'.DS.vel.r')
                    t=read(staname[k]+".subfault"+num+'.DS.vel.t')
                    z=read(staname[k]+".subfault"+num+'.DS.vel.z')
                else: #No! This is the impulse response!
                    r=read(staname[k]+".subfault"+num+'.DS.vel.ri')
                    t=read(staname[k]+".subfault"+num+'.DS.vel.ti')
                    z=read(staname[k]+".subfault"+num+'.DS.vel.zi')
                ntemp,etemp=rt2ne(r[0].data,t[0].data,az[k])
                n=r.copy()
                n[0].data=ntemp/100
                e=t.copy()
                e[0].data=etemp/100
                z[0].data=z[0].data/100
                n=origin_time(n,time_epi,tb)
                e=origin_time(e,time_epi,tb)
                z=origin_time(z,time_epi,tb)
                n.write(staname[k]+".subfault"+num+'.DS.vel.n',format='SAC')
                e.write(staname[k]+".subfault"+num+'.DS.vel.e',format='SAC')
                z.write(staname[k]+".subfault"+num+'.DS.vel.z',format='SAC')
                silentremove(staname[k]+".subfault"+num+'.DS.vel.r')
                silentremove(staname[k]+".subfault"+num+'.DS.vel.t')
                if impulse==True:
                    silentremove(staname[k]+".subfault"+num+'.DS.vel.ri')
                    silentremove(staname[k]+".subfault"+num+'.DS.vel.ti')
                    silentremove(staname[k]+".subfault"+num+'.DS.vel.zi')
        else: #Compute static synthetics
            if okada==False:  #Layered earth model
                
                #this is because when I first wrote this code it processed each
                #source/station pair independently but now that it's vectorized
                #it's does ALL stations in one fell swoop, given the logic it's 
                #easier to just keep this inside the for loop and use the if to
                #run it jsut the first time for all sites
                if k==0:   
                    
                    #initalize output variables
                    staticsSS = zeros((len(d),4))
                    staticsDS = zeros((len(d),4))
                    write_df=True
                    
                    #read GFs file
                    if insar==True:
                        green_file=green_path+model_name+".static."+strdepth+".sub"+subfault+'.insar' #Output dir
                    else: #GPS
                        green_file=green_path+model_name+".static."+strdepth+".sub"+subfault+'.gps' #Output 

                    statics=loadtxt(green_file) #Load GFs
                    Nsites=len(statics) 
                    
                    if len(statics)<1:
                        print('ERROR: Empty GF file')
                        break
                    
                    #Now get radiation pattern terms, there will be 3 terms
                    #for each direction so 9 terms total. THis comes from funtion
                    #dc_radiat() in radiats.c from fk
                    radiation_pattern_ss = zeros((Nsites,9))
                    radiation_pattern_ds = zeros((Nsites,9))
                    
                    rakeSSrad = deg2rad(rakeSS)
                    rakeDSrad = deg2rad(rakeDS)
                    dip_rad = deg2rad(dip)
                    pseudo_strike = deg2rad(az-strike)
                    
                    #Let's do SS first
                    r = rakeSSrad
                    
                    #trigonometric terms following nomenclature used in radiats.c
                    sstk=sin(pseudo_strike) ; cstk=cos(pseudo_strike)
                    sdip=sin(dip_rad) ; cdip=cos(dip_rad)
                    srak=sin(r) ; crak=cos(r)
                    sstk2=2*sstk*cstk ; cstk2=cstk*cstk-sstk*sstk
                    sdip2=2*sdip*cdip ; cdip2=cdip*cdip-sdip*sdip
                    
                    # terms for up component
                    u_dd = 0.5*srak*sdip2*ones(Nsites)
                    u_ds = -sstk*srak*cdip2+cstk*crak*cdip
                    u_ss = -sstk2*crak*sdip-0.5*cstk2*srak*sdip2
                    
                    #terms for r component
                    r_dd = u_dd.copy()
                    r_ds = u_ds.copy()
                    r_ss = u_ss.copy()
                    
                    #terms for t component
                    t_dd = zeros(Nsites)
                    t_ds = cstk*srak*cdip2+sstk*crak*cdip
                    t_ss = cstk2*crak*sdip-0.5*sstk2*srak*sdip2

                    #assemble in one variable
                    radiation_pattern_ss[:,0] = u_dd
                    radiation_pattern_ss[:,1] = u_ds
                    radiation_pattern_ss[:,2] = u_ss
                    
                    radiation_pattern_ss[:,3] = r_dd
                    radiation_pattern_ss[:,4] = r_ds
                    radiation_pattern_ss[:,5] = r_ss
                    
                    radiation_pattern_ss[:,6] = t_dd
                    radiation_pattern_ss[:,7] = t_ds
                    radiation_pattern_ss[:,8] = t_ss
                    
                    
                    #Now radiation pattern for DS
                    r = rakeDSrad
                    
                    #trigonometric terms following nomenclature used in radiats.c
                    sstk=sin(pseudo_strike) ; cstk=cos(pseudo_strike)
                    sdip=sin(dip_rad) ; cdip=cos(dip_rad)
                    srak=sin(r) ; crak=cos(r)
                    sstk2=2*sstk*cstk ; cstk2=cstk*cstk-sstk*sstk
                    sdip2=2*sdip*cdip ; cdip2=cdip*cdip-sdip*sdip
                    
                    # terms for up component
                    u_dd = 0.5*srak*sdip2*ones(Nsites)
                    u_ds = -sstk*srak*cdip2+cstk*crak*cdip
                    u_ss = -sstk2*crak*sdip-0.5*cstk2*srak*sdip2
                    
                    #terms for r component
                    r_dd = u_dd.copy()
                    r_ds = u_ds.copy()
                    r_ss = u_ss.copy()
                    
                    #terms for t component
                    t_dd = zeros(Nsites)
                    t_ds = cstk*srak*cdip2+sstk*crak*cdip
                    t_ss = cstk2*crak*sdip-0.5*sstk2*srak*sdip2
                    
                    #assemble in one variable
                    radiation_pattern_ds[:,0] = u_dd
                    radiation_pattern_ds[:,1] = u_ds
                    radiation_pattern_ds[:,2] = u_ss
                    
                    radiation_pattern_ds[:,3] = r_dd
                    radiation_pattern_ds[:,4] = r_ds
                    radiation_pattern_ds[:,5] = r_ss
                    
                    radiation_pattern_ds[:,6] = t_dd
                    radiation_pattern_ds[:,7] = t_ds
                    radiation_pattern_ds[:,8] = t_ss
                    
                    
                    #Now define the scalng based on magnitude this is variable
                    #"coef" in the syn.c original source code
                    scale = 10**(1.5*Mw+16.1-20) #definition used in syn.c
                    
                    #Scale radiation patterns accordingly
                    radiation_pattern_ss *= scale
                    radiation_pattern_ds *= scale
                    
                    #Now multiply each GF component by the appropriate SCALED
                    #radiation pattern term and add em up to get the displacements
                    # also /100 to convert  to meters
                    up_ss = radiation_pattern_ss[:,0:3]*statics[:,[1,4,7]] 
                    up_ss = up_ss.sum(axis=1) / 100
                    up_ds = radiation_pattern_ds[:,0:3]*statics[:,[1,4,7]] 
                    up_ds = up_ds.sum(axis=1)    / 100                     

                    radial_ss = radiation_pattern_ss[:,3:6]*statics[:,[2,5,8]] 
                    radial_ss = radial_ss.sum(axis=1) / 100
                    radial_ds = radiation_pattern_ds[:,3:6]*statics[:,[2,5,8]] 
                    radial_ds = radial_ds.sum(axis=1) / 100 
                    
                    tangential_ss = radiation_pattern_ss[:,6:9]*statics[:,[3,6,9]] 
                    tangential_ss = tangential_ss.sum(axis=1) / 100
                    tangential_ds = radiation_pattern_ds[:,6:9]*statics[:,[3,6,9]] 
                    tangential_ds = tangential_ds.sum(axis=1) / 100
                    
                    #rotate to neu
                    n_ss,e_ss=rt2ne(radial_ss,tangential_ss,az)
                    n_ds,e_ds=rt2ne(radial_ds,tangential_ds,az)

                    #put in output variables
                    staticsSS[:,0]=n_ss
                    staticsSS[:,1]=e_ss
                    staticsSS[:,2]=up_ss
                    staticsSS[:,3]=beta*ones(Nsites)
                    
                    staticsDS[:,0]=n_ds
                    staticsDS[:,1]=e_ds
                    staticsDS[:,2]=up_ds
                    staticsDS[:,3]=beta*ones(Nsites)
                    
                else:
                    pass
        

            else: #Okada half space solutions
            #SS
                n,e,u=okada_synthetics(strike,dip,rakeSS,ss_length_in_km,ds_length_in_km,xs,ys,
                    zs,lon_sta[k],lat_sta[k],mu_okada)
                savetxt(staname[k]+'.subfault'+num+'.SS.static.neu',(n,e,u,beta),header='north(m),east(m),up(m),beta(degs)')
                #DS
                n,e,u=okada_synthetics(strike,dip,rakeDS,ss_length_in_km,ds_length_in_km,xs,ys,
                    zs,lon_sta[k],lat_sta[k],mu_okada)
                savetxt(staname[k]+'.subfault'+num+'.DS.static.neu',(n,e,u,beta),header='north(m),east(m),up(m),beta(degs)')

    
    if write_df==True and static ==1: #Note to self: stop using 0,1 and swithc to True/False
        
        #Strike slip
        SSdf = df(data=None, index=None, columns=['staname','n','e','u','beta'])
        SSdf.staname=staname
        SSdf.n=staticsSS[:,0]
        SSdf.e=staticsSS[:,1]
        SSdf.u=staticsSS[:,2]
        SSdf.beta=staticsSS[:,3]
        SSdf.to_csv(green_path+'subfault'+num+'.SS.static.neu',sep='\t',index=False,header=False)
        
        DSdf = df(data=None, index=None, columns=['staname','n','e','u','beta'])     
        DSdf.staname=staname
        DSdf.n=staticsDS[:,0]
        DSdf.e=staticsDS[:,1]
        DSdf.u=staticsDS[:,2]
        DSdf.beta=staticsDS[:,3]
        DSdf.to_csv(green_path+'subfault'+num+'.DS.static.neu',sep='\t',index=False,header=False)

    print('Pool  Process %i done working on subfault %i' \
            % (ppool.pid, source[0]))
                        
    # end of make_synthetics_one_subfault



def make_parallel_synthetics_multip(home,project_name,station_file,
                    fault_name,model_name,integrate,static,tsunami,beta,
                    hot_start,time_epi,ncpus,custom_stf,impulse=False,
                    insar=False,okada=False,mu_okada=45e9,NFFT=1024):
    '''
    [New version]
    This routine will take the impulse response (GFs) and pass it into the routine that will
    convovle them with the source time function according to each subfaults strike and dip.
    The result fo this computation is a time series dubbed a "synthetic"
    
    IN:
        home: Home directory
        project_name: Name fo the problem
        station_file: File with coordinates of stations
        fault_name: Name of fault file
        model_Name: Name of Earth structure model file
        integrate: =0 if you want output to be velocity, =1 if you want output to be displacement
        static: =0 if computing full waveforms, =1 if computing only the static field
        hot_start: =k if you want to start computations at k-th subfault, =0 to compute all
        coord_type: =0 if problem is in cartesian coordinates, =1 if problem is in lat/lon
        
    OUT:
        Nothing
    '''
    from numpy import arange,savetxt
    import datetime
    from numpy import loadtxt
    import subprocess
    from shlex import split
    from os import environ
    from multiprocessing import Pool, current_process

    
    station_file=home+project_name+'/data/station_info/'+station_file
    fault_file=home+project_name+'/data/model_info/'+fault_name
    #Time for log file
    now=datetime.datetime.now()
    now=now.strftime('%b-%d-%H%M')
    #First read fault model file
    source=loadtxt(fault_file,ndmin=2)
    
    
    #What parameters are we using?
    if 1:
        out='''Running all processes with:
        home = %s
        project_name = %s
        station_file = %s
        model_name = %s
        integrate = %s
        static = %s
        tsunami = %s
        time_epi = %s
        beta = %d
        custom_stf = %s
        impulse = %s
        insar = %s
        okada = %s
        mu = %.2e
        ''' %(home,project_name,station_file,model_name,str(integrate),\
              str(static),str(tsunami),str(time_epi),beta,custom_stf,impulse,\
              insar,okada,mu_okada)
        print(out)
        
    # Make a list of arguments needed for make_synthetics_one_subfault
    # Only the first argument one_source changes as we loop through all the 
    # subfaults. (Which will be done in parallel with multiprocessing.Pool.)
    
    ksources = list(range(hot_start,source.shape[0]))
    print('Will make synthetics for %i subfaults...' % len(ksources))
    
    all_args = []
    for ksource in ksources:
        one_source = source[ksource,:]  # one row from list of subfaults
        args = (one_source,home,project_name,station_file,\
                            fault_name,model_name,integrate,static,tsunami,\
                            beta, hot_start,time_epi,ncpus,custom_stf,impulse,\
                            insar,okada,mu_okada)
        all_args.append(args)
        
    # now all_args is a list and each element of the list is a tuple with all
    # the arguments needed for make_green_one_subfault.  
    # Farm these out to the requested number of threads, ncpus:
    
    with Pool(processes=ncpus) as pool:
        pool.starmap(make_synthetics_one_subfault, all_args)
    

    # NEW code to store synthetics for all subfaults in a single file
    # for each station.  (Then waveforms for this station can be made using
    # this one file plus the slips for any desired rupture.) 
                

    from obspy import read
    import numpy
    import os

    staname = numpy.genfromtxt(station_file,dtype="U",usecols=0)
    staname = numpy.array(staname,ndmin=1)

    syn_path = home+project_name+'/GFs/synthetics'
    os.system('mkdir -p %s' % syn_path)


    for sta in staname:
        all_synthetics = numpy.empty((len(ksources),6,NFFT), dtype=numpy.float32)
        for ksource in ksources:
            num = str(ksource+1).zfill(4)
            strdepth = '%.4f' % source[ksource,3]
            green_path=home+project_name+'/GFs/dynamic/'+ \
                        "%s_%s.sub%s/" % (model_name,strdepth,num)
            print('+++ Reading from green_path = ',green_path)
            col = 0
            for D in ['SS','DS']:
                for c in ['n','e','z']:
                    fname = green_path + \
                            "/%s.subfault%s.%s.disp.%s" % (sta,num,D,c)
                    synth = read(fname)
                    all_synthetics[ksource,col,:] = synth[0].data
                    col += 1

        fname = syn_path+'/%s.npy' % sta
        numpy.save(fname, all_synthetics)
        print('Saved all synthetics for station %s in %s' % (sta,fname))
        dirnames = home+project_name+'/GFs/dynamic/%s*/%s.subfault*.disp.*' \
                    % (model_name,sta)
        print('Can now delete ', dirnames)
    
         
#Compute GFs for the ivenrse problem            
def inversionGFs(home,project_name,GF_list,tgf_file,fault_name,model_name,
        dt,tsun_dt,NFFT,tsunNFFT,green_flag,synth_flag,dk,pmin,
        pmax,kmax,beta,time_epi,hot_start,ncpus,custom_stf,impulse=False):
    '''
    This routine will read a .gflist file and compute the required GF type for each station
    '''
    from numpy import genfromtxt,where,loadtxt,shape,floor
    from os import remove
    from gc import collect
    from numpy import array
    
    #Read in GFlist and decide what to compute
    gf_file=home+project_name+'/data/station_info/'+GF_list
    stations=genfromtxt(gf_file,usecols=0,skip_header=1,dtype='U')
    stations = array(stations, ndmin=1) # in case only one station
    GF=genfromtxt(gf_file,usecols=[1,2,3,4,5,6,7],skip_header=1,dtype='f8')
    GF = array(GF, ndmin=2) # in case only one station
    fault_file=home+project_name+'/data/model_info/'+fault_name  
    source=loadtxt(fault_file,ndmin=2)
    num_faults=shape(source)[0]
    
    if num_faults/ncpus < 1:
        ncpus=num_faults
        print('Cutting back to ' + str(ncpus) + ' cpus for ' + str(num_faults) + ' subfaults')
            
    # GFs can be computed all at the same time
    station_file='temp.sta'
    try:
        remove(home+project_name+'/data/station_info/'+station_file) #Cleanup
    except:
        pass
    if green_flag==1:
        #decide what GF computation is required for this station
        i=where(GF[:,2]==1)[0]
        if len(i)>0: #Static offset
            print('Static GFs requested...')
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)): #Write temp .sta file
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            static=1
            tsunami=False
            insar=False
            make_parallel_green_multip(home,project_name,station_file,
                        fault_name,model_name,dt,NFFT,static,tsunami,
                        hot_start,dk,pmin,pmax,kmax,ncpus,insar)
        i=where(GF[:,3]==1)[0]
        if len(i)>0 : #displ waveform
            print('Displacememnt GFs requested...')
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)): #Write temp .sta file
                out=str(stations[i[k]])+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            static=0
            tsunami=False
            make_parallel_green_multip(home,project_name,station_file,
                        fault_name,model_name,dt,NFFT,static,tsunami,
                        hot_start,dk,pmin,pmax,kmax,ncpus)
        i=where(GF[:,4]==1)[0]
        if len(i)>0 : #vel waveform
            print('Velocity GFs requested...')
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)): #Write temp .sta file
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            static=0
            tsunami=False
            make_parallel_green_multip(home,project_name,station_file,
                        fault_name,model_name,dt,NFFT,static,tsunami,
                        hot_start,dk,pmin,pmax,kmax,ncpus)
        if tgf_file!=None: #Tsunami
            print('Seafloor displacement GFs requested...')
#            static=0
            static=1
            tsunami=True
            station_file=tgf_file
            make_parallel_green_multip(home,project_name,station_file,
                        fault_name,model_name,dt,NFFT,static,tsunami,
                        hot_start,dk,pmin,pmax,kmax,ncpus)
        i=where(GF[:,6]==1)[0]
        if len(i)>0: #InSAR LOS
            print('InSAR GFs requested...')
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)): #Write temp .sta file
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            static=1
            tsunami=False
            insar=True
            make_parallel_green_multip(home,project_name,station_file,
                        fault_name,model_name,dt,NFFT,static,tsunami,
                        hot_start,dk,pmin,pmax,kmax,ncpus,insar)
            collect()   
    #Synthetics are computed  one station at a time
    if synth_flag==1:
        #Paralell processing
        station_file='temp.sta'
        #Decide which synthetics are required
        i=where(GF[:,2]==1)[0]
        if len(i)>0: #Static offset
            print('Static synthetics requested')
            #Make dummy station file
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)):
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            integrate=0
            static=1
            tsunami=False
            insar=False
            make_parallel_synthetics_multip(home,project_name,station_file,
                fault_name,model_name,integrate,static,tsunami,beta,hot_start,
                time_epi,ncpus,custom_stf,impulse,insar,NFFT=NFFT)
        #Decide which synthetics are required
        i=where(GF[:,3]==1)[0]
        if len(i)>0: #dispalcement waveform
            print('Displacement synthetics requested')
            #Make dummy station file
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)):
                out=str(stations[i[k]])+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            integrate=1
            static=0
            if tgf_file==None: # I am computing for stations on land
                tsunami=False
            else: #I am computing seafloor deformation for tsunami GFs, eventaully
                tsunami=True
            make_parallel_synthetics_multip(home,project_name,station_file,
                fault_name,model_name,integrate,static,tsunami,beta,hot_start,
                time_epi,ncpus,custom_stf,impulse,NFFT=NFFT)
        #Decide which synthetics are required
        i=where(GF[:,4]==1)[0]
        if len(i)>0: #velocity waveform
            print('Velocity synthetics requested')
            #Make dummy station file
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)):
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            integrate=0
            static=0
            if tgf_file==None: # I am computing for stations on land
                tsunami=False
            else: #I am computing seafloor deformation for tsunami GFs, eventaully
                tsunami=True
            make_parallel_synthetics_multip(home,project_name,station_file,
                fault_name,model_name,integrate,static,tsunami,beta,hot_start,
                time_epi,ncpus,custom_stf,impulse,NFFT=NFFT)
        #Decide which synthetics are required
        i=where(GF[:,5]==1)[0]
        if len(i)>0: #tsunami waveform
            print('Tsunami synthetics requested')
            #Make dummy station file
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)):
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            integrate=1
            static=1
            tsunami=True
            station_file=tgf_file
            make_parallel_synthetics_multip(home,project_name,station_file,
                fault_name,model_name,integrate,static,tsunami,beta,hot_start,
                time_epi,ncpus,custom_stf,impulse,NFFT=NFFT)
        #Decide which synthetics are required
        i=where(GF[:,6]==1)[0]
        if len(i)>0: # InSAR LOS
            print('InSAR synthetics requested')
            #Make dummy station file
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)):
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            integrate=0
            static=1
            tsunami=False
            insar=True
            make_parallel_synthetics_multip(home,project_name,station_file,
                fault_name,model_name,integrate,static,tsunami,beta,hot_start,
                time_epi,ncpus,custom_stf,impulse,insar,NFFT=NFFT)
    
                    
  


#Compute GFs for the ivenrse problem            
def teleseismicGFs(home,project_name,GF_list_teleseismic,fault_name,model_name,teleseismic_vel_mod,time_epi,endtime,ncpus):
    '''
    This routine will read a .gflist file and compute the required GF type for each station
    '''
    from numpy import genfromtxt,shape,floor
    from os import remove
    
    #Read in GFlist and decide what to compute
    gf_file=home+project_name+'/data/station_info/'+GF_list_teleseismic
    stations=genfromtxt(gf_file,usecols=0,dtype='U')
    lonlat=genfromtxt(gf_file,usecols=[1,2,])
    fault_file=home+project_name+'/data/model_info/'+fault_name  
    source=genfromtxt(fault_file)
    num_faults=shape(source)[0]
    
    if num_faults/ncpus < 2:
        ncpus=int(floor(num_faults/2.))
        print('Cutting back to ' + str(ncpus) + ' cpus for ' + str(num_faults) + ' subfaults')
    # GFs can be computed all at the same time
    station_file='temp.sta'
    try:
        remove(home+project_name+'/data/station_info/'+station_file) #Cleanup
    except:
        pass
    
    print('Teleseismic GFs requested...')
    f=open(home+project_name+'/data/station_info/'+station_file,'w')
    for k in range(len(stations)): #Write temp .sta file
        out=stations[k]+'\t'+repr(lonlat[k,0])+'\t'+repr(lonlat[k,1])+'\n'
        f.write(out)
    f.close()

    make_parallel_teleseismics_green(home,project_name,station_file,fault_name,model_name,teleseismic_vel_mod,time_epi,endtime,ncpus)
 




                              
                                                        
def run_inversion(home,project_name,run_name,fault_name,model_name,GF_list,G_from_file,G_name,epicenter,
                rupture_speed,num_windows,reg_spatial,reg_temporal,nfaults,beta,decimate,bandpass,
                solver,bounds,weight=False,Ltype=2,target_moment=None,data_vector=None,weights_file=None,
                onset_file=None,GOD_inversion=False):
    '''
    Assemble G and d, determine smoothing and run the inversion
    '''
    from mudpy import inverse as inv
    from mudpy.forward import get_mu_and_area
    from numpy import zeros,dot,array,squeeze,expand_dims,empty,tile,eye,ones,arange,load,size,genfromtxt
    from numpy import where,sort,r_,diag
    from numpy.linalg import lstsq
    from scipy.linalg import norm
    from scipy.sparse import csr_matrix as sparse
    from scipy.optimize import nnls
    from datetime import datetime
    import gc
    from matplotlib import path
    
    

    t1=datetime.now()
    #Get data vector
    if data_vector==None:
        d=inv.getdata(home,project_name,GF_list,decimate,bandpass=bandpass)
    else:
        d=load(data_vector) 
    #Get GFs
    G=inv.getG(home,project_name,fault_name,model_name,GF_list,G_from_file,G_name,epicenter,
                rupture_speed,num_windows,decimate,bandpass,onset_file=onset_file)
    print(G.shape)                
    gc.collect()
    #Get data weights
    if weight==True:
        print('Applying data weights')
        if weights_file==None:
            w=inv.get_data_weights(home,project_name,GF_list,d,decimate)
        else:  # Remember weights are "uncertainties" alrger value is less trustworthy
            w = genfromtxt(weights_file)
            w = 1/w
        
        #apply data weights
        wd=w*d.squeeze()
        
        #get norm after applying weights and normalize weghted data vector
        data_norm = norm(wd)
        wd = wd/data_norm
        
        #reshape for inversion
        wd=expand_dims(wd,axis=1)
        
        #Apply weights to left hand side of the equation (the GFs)
        W=empty(G.shape)
        #W=tile(w,(G.shape[1],1)).T #why this and not a diagonal matrix??? ANSWER: they have the same effect, don;t rememebr why I chose to do it this way
        W=diag(w)
        
        #Normalization effect
        W = W/data_norm
        WG=W.dot(G)

        #Clear up extraneous variables
        W=None
        w=None
        #Define inversion quantities
        x=WG.transpose().dot(wd)
        print('Computing G\'G')
        K=(WG.T).dot(WG)
    else:
        #Define inversion quantities if no weighted
        x=G.transpose().dot(d)
        print('Computing G\'G')
        K=(G.T).dot(G)
    #Get regularization matrices (set to 0 matrix if not needed)
    static=False #Is it jsut a static inversion?
    if size(reg_spatial)>1:
        if Ltype==2: #Laplacian smoothing
            Ls=inv.getLs(home,project_name,fault_name,nfaults,num_windows,bounds)
        elif Ltype==0: #Tikhonov smoothing
            N=nfaults[0]*nfaults[1]*num_windows*2 #Get total no. of model parameters
            Ls=eye(N) 
        elif Ltype==3:  #moment regularization
            N=nfaults[0]*nfaults[1]*num_windows*2 #Get total no. of model parameters
            Ls=ones((1,N))
            #Add rigidity and subfault area
            mu,area=get_mu_and_area(home,project_name,fault_name,model_name)
            istrike=arange(0,N,2)
            Ls[0,istrike]=area*mu
            idip=arange(1,N,2)
            Ls[0,idip]=area*mu
            #Modify inversion quantities
            x=x+Ls.T.dot(target_moment)
        else:
            print('ERROR: Unrecognized regularization type requested')
            return
        Ninversion=len(reg_spatial)
    else:
        Ls=zeros(K.shape)
        reg_spatial=array([0.])
        Ninversion=1
    if size(reg_temporal)>1:
        Lt=inv.getLt(home,project_name,fault_name,num_windows)
        Ninversion=len(reg_temporal)*Ninversion
    else:
        Lt=zeros(K.shape)
        reg_temporal=array([0.])
        static=True
    #Make L's sparse
    Ls=sparse(Ls)
    Lt=sparse(Lt)
    #Get regularization tranposes for ABIC
    LsLs=Ls.transpose().dot(Ls)
    LtLt=Lt.transpose().dot(Lt)
    #inflate
    Ls=Ls.todense()
    Lt=Lt.todense()
    LsLs=LsLs.todense()
    LtLt=LtLt.todense()
    #off we go
    dt=datetime.now()-t1
    print('Preprocessing wall time was '+str(dt))
    print('\n--- RUNNING INVERSIONS ---\n')
    ttotal=datetime.now()
    kout=0
    for kt in range(len(reg_temporal)):
        for ks in range(len(reg_spatial)):
            t1=datetime.now()
            lambda_spatial=reg_spatial[ks]
            lambda_temporal=reg_temporal[kt]
            print('Running inversion '+str(kout+1)+' of '+str(Ninversion)+' at regularization levels: ls ='+repr(lambda_spatial)+' , lt = '+repr(lambda_temporal))
            if static==True: #Only statics inversion no Lt matrix
                Kinv=K+(lambda_spatial**2)*LsLs
                Lt=eye(len(K))
                LtLt=Lt.T.dot(Lt)
            else: #Mixed inversion
                Kinv=K+(lambda_spatial**2)*LsLs+(lambda_temporal**2)*LtLt
            if solver.lower()=='lstsq':
                sol,res,rank,s=lstsq(Kinv,x)
            elif solver.lower()=='nnls':
                x=squeeze(x.T)
                try:
                    sol,res=nnls(Kinv,x)
                except:
                    print('+++ WARNING: No solution found, writting zeros.')
                    sol=zeros(G.shape[1])
                x=expand_dims(x,axis=1)
                sol=expand_dims(sol,axis=1)
            else:
                print('ERROR: Unrecognized solver \''+solver+'\'')

            #Force faults outside a polygon to be zero
#            print('WARNING: Using fault polygon to force solutions to zero')
#            #load faulkt
#            fault=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
#            polygon=genfromtxt('/Users/dmelgarm/Oaxaca2020/etc/zero_fault.txt')
#            polygon=path.Path(polygon)
#            i=where(polygon.contains_points(fault[:,1:3])==False)[0]
#            i=sort(r_[i*2,i*2+1])
#            N=nfaults[0]*2
#            i=r_[i,i+N,i+2*N,i+3*N]
#            sol[i]=0
            
            #Compute synthetics
            ds=dot(G,sol)
            
            #Get stats
            L2,Lmodel=inv.get_stats(Kinv,sol,x,Ls)
            VR,L2data=inv.get_VR(home,project_name,GF_list,sol,d,ds,decimate,WG,wd)
            #VR=inv.get_VR(WG,sol,wd)
            #ABIC=inv.get_ABIC(WG,K,sol,wd,lambda_spatial,lambda_temporal,Ls,LsLs,Lt,LtLt)
            ABIC=inv.get_ABIC(G,K,sol,d,lambda_spatial,lambda_temporal,Ls,LsLs,Lt,LtLt)
            #Get moment
            Mo,Mw=inv.get_moment(home,project_name,fault_name,model_name,sol)
            #If a rotational offset was applied then reverse it for output to file
            if beta !=0:
                sol=inv.rot2ds(sol,beta)
            #Write log
            inv.write_log(home,project_name,run_name,kout,rupture_speed,num_windows,lambda_spatial,lambda_temporal,beta,
                L2,Lmodel,VR,ABIC,Mo,Mw,model_name,fault_name,G_name,GF_list,solver,L2data)
            #Write output to file
            if GOD_inversion==True:
                num=str(kout).rjust(4,'0')
                np.save(home+project_name+'/output/inverse_models/'+run_name+'.'+num+'.syn.npy',ds)
                inv.write_synthetics_GOD(home,project_name,run_name,GF_list,ds,kout,decimate)
            else:
                inv.write_synthetics(home,project_name,run_name,GF_list,G,sol,ds,kout,decimate)
            inv.write_model(home,project_name,run_name,fault_name,model_name,rupture_speed,num_windows,epicenter,sol,kout,onset_file=onset_file)
            kout+=1
            dt1=datetime.now()-t1
            dt2=datetime.now()-ttotal
            print('... inversion wall time was '+str(dt1)+', total wall time elapsed is '+str(dt2))
