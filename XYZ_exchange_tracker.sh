
#!/bin/sh
SCRIPT="XYZ_user_var"
F90="$SCRIPT"".f90"
EXE="$SCRIPT"".exe"


XYZ=$1
natom=`head -1 $XYZ`
SYSTEM=`echo ${XYZ/".xyz"/""}`
CP2K_IN='./md.in'
if [ -f "$CP2K_IN" ]; then
    ABC=`grep ABC md.in | head -1 | awk '{print $3, $4, $5}'`
else 
    echo "Enter ORTHORHOMBIC cell parameters: "
    read ABC
fi

gfortran -ofast -fcheck=bound $F90 -o $EXE

./$EXE $ABC $natom $SYSTEM $XYZ

echo ""
echo ""

