import h5py

def read_data(filename,dataname,  slice = 'none'):

    with h5py.File(filename, "r") as f:

        nslice=int(slice[7:])

        # print(nslice,type(nslice))



        if slice[0:6]=='kslice':
            data_read=f[dataname][nslice,:,:]
        elif slice[0:6] =='jslice':
            data_read=f[dataname][:,nslice,:]
        elif slice[0:6]=='islice':
            data_read=f[dataname][:,:,nslice]
        else:
            data_read=f[dataname][:]

        f.close()

        # print(slice[0:6])
        # 

        # print(filename,dataname)
        # print(data_read)
        return data_read

        # data_list=f.keys()
        # print(data_list,type(data_list))
    # Print all root level object names (aka keys) 
    # these can be group or dataset names 
    # print("Keys: %s" % f.keys())
    # get first object name/key; may or may NOT be a group
    # a_group_key = list(f.keys())[0]

    # # get the object type for a_group_key: usually group or dataset
    # print(type(f[a_group_key])) 

    # # If a_group_key is a group name, 
    # # this gets the object names in the group and returns as a list
    # data = list(f[a_group_key])

    # # If a_group_key is a dataset name, 
    # # this gets the dataset values and returns as a list
    # data = list(f[a_group_key])
    # # preferred methods to get dataset values:
    # ds_obj = f[a_group_key]      # returns as a h5py dataset object
    # ds_arr = f[a_group_key][()]  # returns as a numpy array

    # print(ds_obj)
    # print(ds_arr)