find . -type f -exec md5sum {} + | sort | awk 'BEGIN{last=""} {if($1==last) {print $2} else {last=$1}}' | xargs -d '\n' rm
find . -type d -empty -delete

