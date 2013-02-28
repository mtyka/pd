rm `find -name "pd*.py"` -f
rm `find -name "*.py" | grep -v run` -f
rm `find -name "*so"` -f
#rm `find -name ".svn*"` -rf
rm `find -name "inout*"` -rf
rm `find -name "*py[co]"` -rf
rm `find -name "output.*"` -f
rm `find -name "*.la"` -f
rm `find -name "*.0"` -f
rm *.log -f
