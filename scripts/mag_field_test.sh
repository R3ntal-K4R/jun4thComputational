gcc magnetic_field_test.c -o magnetic_field_test -lm
if [ $? -eq 0 ]; then
echo "Compiled succesfully"
./magnetic_field_test
./mag_field_test_plot.sh
fi
if [ $? -ne 0 ]; then
echo "Did not compile"
fi
