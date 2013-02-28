echo Performing python modifications ...
python ..\src\tools\modifypd.py  ../src/pd/pd.py
echo Copying executables to bin folder...
copy ..\src\pd\pd.py ..\bin
IF EXIST ..\bin\_pd.dll copy ..\bin\_pd.dll ..\bin\_pd.pyd
IF EXIST ..\bin\_arcus.dll copy ..\bin\_arcus.dll ..\bin\_arcus.pyd
IF EXIST ..\bin\_example.dll copy ..\bin\_example.dll ..\bin\_example.pyd