# leechem: Python code for the Leech SBEM Datasets - ReadMe

## 1. Introduction
The *leechem* repository provides python code to retrieve data from SBEM (serial-block face electron microscopy) database.  

- The **Database**: Contains SBEM reconstruction of the canonical ganglion 10 from the leech. The focus has in particular been on the reconstruction of presynaptic projections to neuron DE3 (tree id 444 in the database), the dorsal excitatory motorneuron 3.

- The **Python code**: Accesses the data collected during the reconstruction and provides for instance x,y,z coordinates on the different measurement points (nodes) which are annotated with their anatomical descriptors (e.g. "soma" or "presynaptic terminal"). The code accesses the data, but also provides functions to perform data analysis. The most important methods are to be found in the sbemdb class within the sbemdb.py python file.

## 2. Installation

### Python Packages

To run the python code, the relevant python packages need to be installed. See for instance the following lines in sbemdb.py:

```
import numpy as np
import sqlite3
import errno
import os
```
This modules need to be installed before running the code, and can be most conveniently done using **pip** from the terminal.
For an introduction to pip, [see here.](https://pip.pypa.io/en/stable/installing/)

## 3. Code

### Python code - A very brief overview

The python code relies mainly on sqlite to retrieve the data from the database. In an SQL environment, the typical process is
to create a database object and then a so called cursor to retrieve with sql specific clauses the desired information.
For instance:

```
db = SBEMDB() #creating a database object
c = db.db.cursor() #creating the cursor c
c.execute("select name from sqlite_master where type = 'table';") #using the cursor c with a clause to retrieve all table names
table_names = c.fecthall() #retrieve all results from c
print(table_names) #print results
```
The code in sbemdb.py has been written, however, in a way that sql specific calls are minimized. 

The function *nodexyz* for instance:

```
def nodexyz(self, where):
        '''NODEXYZ - Get a map of node positions
        (x, y, z, nid) = NODEXYZ(whereclause) performs a SELECT on the NODES 
        with the given where clause. For instance,
        nnn = db.nodes('tid==444 and typ==6') finds
        all the postsynaptic terminals on tree #444.
        Results X, Y, Z are expressed in um and are corrected for the gapshift.
        ID is node ID'''
        c = self.db.cursor()
        c.execute('select nid, x, y, z from nodes where %s' % where)
        rows = c.fetchall()
       (...)
        return (x, y, z, nid)
```

allows the user to retrieve the 3D coordinates of a note and its id ("nid") via for instance "nodexyz('tid==444')", while
the object oriented code creates a cursor and accesses the database hidden from the user.

The most convenient way to get started with the code might be by way of the "demo.py" program.



