
for i in $*; do echo $i; done|cut -d '-' -f 1|sort -u |
while read line; do 
	cat $line*.dump > $line.all
done


