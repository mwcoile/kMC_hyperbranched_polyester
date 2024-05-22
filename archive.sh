# archives current version of code as a previous version that can be easily restored or viewed

CurrentDate=$(date +"%y%m%d")
Comment="backupaftersplitting"

files=("homo.cpp" "input.txt" "Makefile" "inputs.cpp" "inputs.h" "chainUpdate.cpp" "chainUpdate.h" "molecularWeight.cpp" "molecularWeight.h" "analysis.cpp" "analysis.h" "tests.cpp" "tests.h")
mkdir "./OldVersions/${CurrentDate}$Comment"
path="./OldVersions/${CurrentDate}$Comment"

for i in "${files[@]}"
do
    echo "$i"
    cp $i "${path}/$i"
done

