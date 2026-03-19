ls `find -name '*.sh'` `find -name '*.py'` `find -name '*.md'` `find -name '*.json'` | egrep -v 'boltz|ses_' > list_scripts.txt
