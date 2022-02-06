# pysrcor 

Python3 version of the ["srcor" function](https://idlastro.gsfc.nasa.gov/ftp/pro/idlphot/srcor.pro) in [NASA's IDL Astronomy User's Library](https://idlastro.gsfc.nasa.gov/contents.html)

A quick and easy function to correlate source positions from two catalogs.

## Install 
The easiest way is using pip:

>pip install pysrcor

## Usage
>from pysrcor.pysrcor import Cat

>ra1, dec1 = [145.4354343, 150.234245], [-27.23423, -30.324233]

>ra2, dec2 = [0.003423, 145.4355343, 150.234235], [10.32432, -27.23423, -30.324233]

>cat = Cat(ra1, dec1, ra2, dec2)

>id1, id2, dis = cat.match(rad=1, opt=2, silent=False)

id1 and id2 are the matched source indexes (starting from 0) in catalogs 1 & 2
dis is the distances of matched sources in arcsec

## Contact
Post any comments/suggestions/questions under "Issues" 
or email to
gyang206265@gmail.com (Guang Yang)
