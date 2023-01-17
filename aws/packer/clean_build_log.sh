uid=$(date +%s)
grep -v amazon-ebs\:....$ logs/build.log.tmp > logs/build.${uid}.log
rm build.log.tmp

