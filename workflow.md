KinMut2
======


To install this workflow you will need a working rbbt installation. Please find
documentation on getting a working rbbt 
instalation [here](http://mikisvaz.github.io/rbbt/tutorial/install/).

In short it involves configuring the base system with some packages needed
by ruby, installing the necessary ruby gems (in particular the rbbt framework)
and then bootstrapping the system (see the instructions for more details on
file_servers)

An alternative version is using the [rbbt-image
package](https://github.com/mikisvaz/rbbt-image) that can construct a provision
script to be used for local instalations as well as producing docker or vagrant
images. 

The simplest way to use it is calling the `predict` task of the workflow, which
limits itself to the KinMut predictions. The predict_all method will call the
StructurePPi and DbNSFP workflow that have a very lengthy boostrapping process,
so its recomended to configure these as `remote workflows` if one is interested
in this option.

Note that while the code for the web server is also provided, part of it
depends on having access to the postgress databases that have not been included
in this repository and will not render any information. 

