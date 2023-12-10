# convert all dos file under current directory to unix format

files=`find . -type f -exec file {} \; | grep CRLF | cut -d: -f1`

for file in $files; do
    dos2unix $file
done
