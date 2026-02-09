#! /bin/bash 
echo
echo --------------------------------
echo writing files for simulation ...
echo --------------------------------

# cp size of the mechanism first
echo ... cp size.dum in $1.mech
cp size.dum $1.mech

# assume gas phase species always exist
echo ... cat gasspe.dum in $1.mech
cat $1.mech gasspe.dum > temp.mech
mv temp.mech $1.mech

# add particle phase species, if any
nlin=$(cat partspe.dum | wc -l)
echo number of species in part. phase: $nlin
if [ $nlin -gt 1 ];
then
echo ... cat partspe.dum '>' $1.mech
cat $1.mech partspe.dum > temp.mech
mv temp.mech $1.mech
fi

# add wall phase species, if any
nlin=$(cat wallspe.dum | wc -l)
echo number of species in wall phase: $nlin
if [ $nlin -gt 1 ];
then
echo ... cat wallspe.dum '>' $1.mech
cat $1.mech wallspe.dum > temp.mech
mv temp.mech $1.mech
fi

# add reactions
echo ... cat reactions.dum '>' $1.mech
cat $1.mech reactions.dum > temp.mech
mv temp.mech $1.mech
