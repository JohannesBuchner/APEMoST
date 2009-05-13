
if [ "$CATCOMMAND" == "" ]; then
	cat="cat"
else
	cat="$CATCOMMAND"
fi
for i in $*; do echo $i; done|cut -d '-' -f 1|sort -u |
while read line; do 
	paste -d '\n' $line*.dump|$cat > $line.all
done

