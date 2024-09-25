This code is to create initial data sets for a boson gas i.e. a complex scalar field.
First we create data with random values (either cuted pseudo-gaussian or thermal), then 
fourier-trasnform them to physical space and cut them to a ball insida a box. 
We then create a hdf5 file to import to the code for evolution.
