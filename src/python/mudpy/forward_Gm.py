'''
D. Melgar 02/2014

Forward modeling routines

Modified version to accumulate G*m
@rjleveque 2/2022

'''

# Things we need that haven't changed from forward.py:
from mudpy.forward import load_fakequakes_synthetics, read_fakequakes_hypo_time
from mudpy.forward import tshift_trace, build_source_time_function
from mudpy.forward import write_fakequakes_waveforms, gauss_reconvolution


def waveforms_fakequakes_Gm(home,project_name,fault_name,rupture_list,GF_list,
                model_name,run_name,dt,NFFT,
                source_time_function='dreger',zeta=0.2,
                stf_falloff_rate=4.0,rupture_name=None,
                epicenter=None,time_epi=None, hot_start=0):
    '''
    To supplant waveforms_matrix() it needs to include resmapling and all that jazz...
    
    Instead of doing matrix multiplication one stations at a time, do it for all stations
    
    This routine will take synthetics and apply a slip dsitribution. It will delay each 
    subfault by the appropriate rupture time and linearly superimpose all of them. Output
    will be one sac waveform file per direction of motion (NEU) for each station defined in the
    station_file. Depending on the specified rake angle at each subfault the code will compute 
    the contribution to dip and strike slip directions. It will also compute the moment at that
    subfault and scale it according to the unit amount of momeent (1e15 N-m)
    
    IN:
        home: Home directory
        project_name: Name of the problem
        rupture_name: Name of rupture description file
        station_file: File with coordinates of stations
        model_Name: Name of Earth structure model file
        integrate: =0 if you want output to be velocity, =1 if you want output to de displacement
       
    OUT:
        Nothing
    '''
    from numpy import genfromtxt,array
    import datetime
    import time
    
    print('Solving for kinematic problem(s)')
    #Time for log file
    now=datetime.datetime.now()
    now=now.strftime('%b-%d-%H%M')
    
    #load source names
    if rupture_list==None:
        #all_sources=array([home+project_name+'/forward_models/'+rupture_name])   
        all_sources=array([rupture_name])  
    else:
        all_sources=genfromtxt(home+project_name+'/data/'+rupture_list,dtype='U')
        all_sources = array(all_sources, ndmin=1)  # in case only 1 entry
    
    #Load all synthetics
    print('... loading all synthetics into memory')
    Nss,Ess,Zss,Nds,Eds,Zds=load_fakequakes_synthetics(home,project_name,
                            fault_name,model_name,GF_list,
                            G_from_file=False,G_name=project_name)
    print('... ... done')
    
    # Need to know how many sites
    station_file=home+project_name+'/data/station_info/'+GF_list
    staname=genfromtxt(station_file,dtype="U",usecols=0)
    staname = array(staname, ndmin=1) # in case only one station
    Nsta=len(staname)
    
    #Now get the impulse response G for all sites and all subfaults
    #Gimpulse_all = Nss,Ess,Zss,Nds,Eds,Zds=mseed2matrix(Nss,,Ess,Zss,Nds,Eds,Zds)
    
    
    #Now loop over rupture models
    for ksource in range(hot_start,len(all_sources)):
        print('... solving for source '+str(ksource)+' of '+str(len(all_sources)))
        rupture_name=all_sources[ksource]
        print('... ' + rupture_name)
        
        if rupture_list!=None:
            #Get epicentral time
            epicenter,time_epi=read_fakequakes_hypo_time(home,project_name,rupture_name)
            forward=False
        else:
            forward=True #This controls where we look for the rupture file
        
        # Compute G*m and return as waveforms:

        t1=time.time()
        
        waveforms = get_fakequakes_G_times_m (Nss,Ess,Zss,Nds,Eds,Zds,home,
                        project_name,rupture_name,time_epi,GF_list,
                        epicenter,NFFT,source_time_function,
                        stf_falloff_rate,zeta,forward=forward)
        t2=time.time()
                        
        print('... ... slip rate convolutions completed in {:.1f}s'.format(t2-t1))
        
        #Write output
        write_fakequakes_waveforms(home,project_name,rupture_name,waveforms,GF_list,NFFT,time_epi,dt)
        
    
def stream2matrix(Nss,Ess,Zss,Nds,Eds,Zds,Nsta):
    '''
    
    Converts stream objects to a properly formatted G matrix

    Parameters
    ----------
    Nss,Ess ... : Stream objects with synthetics for each component of motion and
                    for the ss and ds rake angles


    Nsta : int, number of stations being processed


    Returns
    -------
    G = [ Nss_sta1_sf1 Nds sta1_ sf1     Nss_sta1_sf2 Nds sta1_ sf2 ...
          Ess_sta1_sf1 Eds sta1_ sf1     Nss_sta1_sf3 Nds sta1_ sf3 ...
          Zss_sta1_sf1 Zds sta1_ sf1     Nss_sta1_sf3 Nds sta1_ sf3 ...
          Nss_sta2_sf1 Nds sta2_ sf1
          Ess_sta2_sf1 Eds sta2_ sf1
          Zss_sta2_sf1 Zds sta2_ sf1
          ...

    '''    
    
    from numpy import ones
    
    #How many time points ins ytnethics
    Npts = Nss[0].stats.npts
    
    #How many sources?
    Nsources = int(len(Nss) / Nsta)
    
    #Initalize G matrix
    G = ones((Npts*Nsta*3,2*Nsources))
    
    k=0
    row_start = 0
    row_end = Npts
    for ksta in range(Nsta):
        for ksub in range(Nsources):
            
            G[row_start:row_end,2*ksub] = Nss[k].data
            G[row_start:row_end,2*ksub+1] = Nds[k].data
            
            G[row_start+Npts:row_end+Npts,2*ksub] = Ess[k].data
            G[row_start+Npts:row_end+Npts,2*ksub+1] = Eds[k].data
    
            G[row_start+2*Npts:row_end+2*Npts,2*ksub] = Zss[k].data
            G[row_start+2*Npts:row_end+2*Npts,2*ksub+1] = Zds[k].data
            
            k += 1
    
        row_start += 3*Npts
        row_end += 3*Npts
    
    return G
        





def get_fakequakes_G_times_m(Nss,Ess,Zss,Nds,Eds,Zds,home,
                project_name,rupture_name,time_epi,GF_list,epicenter,NFFT,
                source_time_function,stf_falloff_rate,zeta=0.2,
                forward=False,reconvolution=False,old_stf='prem_i_2s'):
    '''
    Directly compute G*m, synthetics shifted in time multiplied by slip vector.
    
    IN:
        home: Home directory location
        project_name: Name of the problem
        fault_name: Name of fault description file
        model_name: Name of velocity structure file
        GF_list: Name of GF control file
        epicenter: Epicenter coordinates
        
        Needs updating
        
    OUT:
        Gm: Product of GF matrix and slip vector m.
    '''
    
    from numpy import genfromtxt,convolve,where,zeros,arange,unique,array
    import gc


    if forward==True:
        source=genfromtxt(home+project_name+'/forward_models/'+rupture_name)
    else:
        source=genfromtxt(home+project_name+'/output/ruptures/'+rupture_name)
    rise_times=source[:,7]
    rupture_onset=source[:,12]
    
    #How many unique faults?
    Nfaults=len(unique(source[:,0]))
    
    #How many subfaults are non-zero?
    i_non_zero=where(rise_times>0)[0]
    N_non_zero=len(i_non_zero)
    
    #Stations
    station_file=home+project_name+'/data/station_info/'+GF_list
    staname=genfromtxt(station_file,dtype="U",usecols=0)
    staname = array(staname, ndmin=1) # in case only one station
    Nsta=len(staname)
    
    #Get slip model vector
    m=zeros((N_non_zero*2,1))
    iss=arange(0,len(m),2)
    ids=arange(1,len(m),2)
    m[iss,0]=source[i_non_zero,8]
    m[ids,0]=source[i_non_zero,9]
    
    #Initalize Gm vector for product G*m
    Gm = zeros(NFFT*3*Nsta)
    
    # Compute two columns of G matrix as in get_G_and_m version, but
    # instead of placing synthetics in correct place in Gmatrix,
    # we will then multiply by corresponding components of m and accumulate.
    
    matrix_pos=0 #tracks where in matrix synths are placed
    read_start=0  #Which trace to start reading from
    for ksta in range(Nsta):
        print('... ... working on station '+str(ksta+1)+' of '+str(Nsta))
        
        for ksource in range(len(i_non_zero)):

            if ksource % 50 == 0:
                print('... ... ... subfault %d of %d' % (ksource,len(i_non_zero)))
            
            #Get synthetics
            nss=Nss[read_start+i_non_zero[ksource]].copy()
            ess=Ess[read_start+i_non_zero[ksource]].copy()
            zss=Zss[read_start+i_non_zero[ksource]].copy()
            nds=Nds[read_start+i_non_zero[ksource]].copy()
            eds=Eds[read_start+i_non_zero[ksource]].copy()
            zds=Zds[read_start+i_non_zero[ksource]].copy()
            #Delay synthetics by rupture onset
            tdelay=rupture_onset[i_non_zero[ksource]]
            nss,ess,zss,nds,eds,zds=tshift_trace(nss,ess,zss,nds,eds,zds,tdelay,time_epi,NFFT)
            #Convolve with source time function
            rise=rise_times[i_non_zero[ksource]]
            #Make sure rise time is a multiple of dt
            dt=nss.stats.delta
            rise=round(rise/dt)*nss.stats.delta
            if rise<(2.*dt): #Otherwise get nan's in STF
                rise=2.*dt
            
            total_time=NFFT*dt-dt
            
            if reconvolution==False: #Do a straight up convolution with new slip rate function
                
                #Build the new slip rate function
                t_stf,stf=build_source_time_function(rise,dt,total_time,stf_type=source_time_function,zeta=zeta,dreger_falloff_rate=stf_falloff_rate)

                nss.data=convolve(nss.data,stf)[0:NFFT]
                ess.data=convolve(ess.data,stf)[0:NFFT]
                zss.data=convolve(zss.data,stf)[0:NFFT]
                nds.data=convolve(nds.data,stf)[0:NFFT]
                eds.data=convolve(eds.data,stf)[0:NFFT]
                zds.data=convolve(zds.data,stf)[0:NFFT]
            else: #reconvovle by using old stf
            

                
                #Need old STF
                old_rise_time = 2.1 #This is the prem_i_2s hard coded value
                time_offset_gauss=100 #Hard coded for now. Long enough for any conceivable rise time
                t_stf,old_stf=build_source_time_function(old_rise_time,dt,total_time,stf_type='gauss_prem_i_2s',time_offset_gauss=time_offset_gauss,scale=True,scale_value=1.0)
                
                #Now the new STF
                t_stf,new_stf=build_source_time_function(rise,dt,total_time,stf_type='gauss_prem_i_2s',time_offset_gauss=time_offset_gauss,scale=True,scale_value=1.0,quiet=True)
                
                #Reconvolutions
                nss.data=gauss_reconvolution(nss.data,new_stf,old_stf)
                ess.data=gauss_reconvolution(ess.data,new_stf,old_stf)
                zss.data=gauss_reconvolution(zss.data,new_stf,old_stf)
                nds.data=gauss_reconvolution(nds.data,new_stf,old_stf)
                eds.data=gauss_reconvolution(eds.data,new_stf,old_stf)
                zds.data=gauss_reconvolution(zds.data,new_stf,old_stf)
                


            # Multiply these columns by corresponding slips m and accumulate
            # in product G*m:
            
            # Note: need to accumulate using both strike and dip components of
            # slip, and the rows are partitioned into 3 parts for 
            # N,E,Z components of waveform
            
            Gm[matrix_pos:matrix_pos+NFFT] += nss.data * m[2*ksource]
            Gm[matrix_pos:matrix_pos+NFFT] += nds.data * m[2*ksource+1]
            Gm[matrix_pos+NFFT:matrix_pos+2*NFFT] += ess.data * m[2*ksource]
            Gm[matrix_pos+NFFT:matrix_pos+2*NFFT] += eds.data * m[2*ksource+1]
            Gm[matrix_pos+2*NFFT:matrix_pos+3*NFFT] += zss.data * m[2*ksource]
            Gm[matrix_pos+2*NFFT:matrix_pos+3*NFFT] += zds.data * m[2*ksource+1]
            
            
            #Having some sort of memory issue delete streams and garbage collect to avoid it
            #del nss,ess,zss,nds,eds,zds
            #gc.collect()
            
        matrix_pos+=3*NFFT
        read_start+=Nfaults
        
    return Gm


