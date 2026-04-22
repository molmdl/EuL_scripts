find /path/to/directory -type f -exec md5sum {} + | sort | awk 'BEGIN{last=""} {if($1==last) {print $2} else {last=$1}}' | xargs -d '\n' rm
