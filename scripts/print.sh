while :
do
   echo Enter Line
   read linenum
   echo "Line:" $linenum
   #sed -ne ${linenum}p $1
   sed -n "${linenum},${p;q;}" $1
   sed -n "${linenum},${p;q;}" $2
   #sed -ne ${linenum}p $2
   #awk "FNR==$linenum {print FILENAME, '$0'}" $1 $2
done
