
awk '$1 == 2 && $2 == 2 && $3 ==2' Bfieldn400.dat > origBline.dat



awk '$1 == 20 && $3 == 20' Bfieldn400.dat > origBline.dat &
awk '$1 == 20 && $3 == 20' mag_out.dat > mag_out_Bline.dat &


plot 'origBline.dat' u 2:4 w l
rep 'mag_out_Bline.dat' u 2:4 w l


plot 'origBline.dat' u 2:5 w l
rep 'mag_out_Bline.dat' u 2:5 w l

plot 'origBline.dat' u 2:6 w l
rep 'mag_out_Bline.dat' u 2:6 w l



# if [ $# -ne 4 ]
# then
# echo Usage: $0 x y z Bfile.dat. Use "*" for scanning variable/
# fi

# i=$1
# j=$2
# k=$3
# Bfile=$4

# #echo awk "$1 == x && $2 == y && $3 ==z" $Bfile
# #awk '$1 == x && $2 == y && $3 ==z' $Bfile

# #gawk -v '$1=$i -v y=$j -v z=$k' $Bfile


# #awk '$1 == 2 && $2 == 2 && $3 ==2' Bfield.da


# x=10
# y=30
# text="Total is : "
# awk -v a=$x -v b=$y -v c="$text" 'BEGIN {ans=a+b; print $1 $2 $3}'