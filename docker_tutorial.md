![MaSIF banner and concept](https://raw.githubusercontent.com/LPDI-EPFL/masif/master/img/Concept-01.png)

# Docker tutorial for MaSIF.

## Installation

```
docker pull pablogainza/masif:latest
docker run -it pablogainza/masif
```
You now start a local container with MaSIF. The first step is to clone the masif\_molecule repository:

```
root@b30c52bcb86f:/masif# git clone https://github.com/jomimc/masif_molecule.git
```

## MaSIF pipeline

### Running the pipeline on a protein from a PDB file

Go into the MaSIF-molecule source directory. 
```
cd masif_molecule/source
```

Running the pipeline on the test data provided "data/2sic.pdb". To extract a single protein (chain) use the "-c" option. Specify output folder with "-o". Run "python main.py --help" to see all options. In this example we compute all possible features and patches with default parameters.

```
python main.py ../data/2sic.pdb -c E -o test_data --hphob --hbond --patches
```

If you want to run a prediction on multiple chains you can run without the "-c" option:

```
python main.py ../data/2sic.pdb -o test_data --hphob --hbond --patches
```


### Running the pipeline on a nanoparticle from a PQR file

Go into the MaSIF-molecule source directory. 
```
cd masif_molecule/source
```

Running the pipeline on the test data provided "data/np.pqr". To extract a single protein (chain) use the "-c" option. Specify output folder with "-o". Run "python main.py --help" to see all options. In this example we compute all possible features and patches with default parameters.

```
python main.py ../data/Np.pdb -o test_data --patches
```

