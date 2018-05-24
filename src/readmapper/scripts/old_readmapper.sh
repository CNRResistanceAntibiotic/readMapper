#!/bin/bash
#Richard Bonnet
#1/07/2017


#usage message
usage() { echo "Usage: $0 [-s sample file] [-r reads directory] [-d work directory] [-b databases directory] [-i initial of the user] <-F Overwrite output directory (Default=False)>" 1>&2; exit 1; }

#get application directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#parse args
force='0'
while getopts "s:r:d:b:i:VF" opt; do
    case "${opt}" in
        s)
            sampleFile=${OPTARG}
            ;;
        r)
            readsDir=${OPTARG}
            ;;
        d)
            wkDir=${OPTARG}
            ;;

	    b)
            dbDir=${OPTARG}
            ;;          
                        
        i)  
            initial=${OPTARG}
            ;;
        F)
            force='1'
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

#case of one or more args lack
if [ -z "${sampleFile}" ] || [ -z "${readsDir}" ] || [ -z "${wkDir}" ] || [ -z "${dbDir}" ] || [ -z "${initial}" ]; then
    usage
fi

#get folders/files path
sampleFile=$(readlink -f ${sampleFile})
readsDir=$(readlink -f ${readsDir})
wkDir=$(readlink -f ${wkDir})
dbDir=$(readlink -f ${dbDir})
setFile=$(readlink -f ${dbDir}/setting.txt)

#print folders/files path
echo -e "\nSample file: ${sampleFile}"
echo    "Reads directory: ${readsDir}"
echo    "Work directory: ${wkDir}"
echo    "DataBase directory: ${dbDir}"
echo    "Setting file: ${setFile}"
echo -e "Application run at : ${DIR}\n"

if [ ! -e "${sampleFile}" ] || [ ! -e "${readsDir}" ] || [ ! -e "${wkDir}" ] || [ ! -e "${dbDir}" ] ; then
    usage
fi

echo "----------------"
echo "PREPARE THE JOBS"
echo "----------------"
#echo $force

if [ "${force}" == '1' ];then
    echo -e "\nForce the preparation : \n"
    echo -e "${DIR}/prepare_mapping.py \n\t -sf ${sampleFile} \n\t -rd ${readsDir} \n\t -wd ${wkDir} \n\t -in ${initial} \n\t -set ${setFile} \n\t -F \n"
	${DIR}/prepare_mapping.py -sf ${sampleFile} -rd ${readsDir} -wd ${wkDir} -in ${initial} -set ${setFile} -F

else
    echo -e "\nRun the preparation : \n"
    echo -e "${DIR}/prepare_mapping.py \n\t -sf ${sampleFile} \n\t -rd ${readsDir} \n\t -wd ${wkDir} \n\t -in ${initial} \n\t -set ${setFile} \n"
	${DIR}/prepare_mapping.py -sf ${sampleFile} -rd ${readsDir} -wd ${wkDir} -in ${initial} -set ${setFile}
fi

echo ""
echo "-----------------"
echo "STARTING THE JOBS"
echo "-----------------"

echo -e "\nRun the manager : \n"
echo -e "${DIR}/manager.py  \n\t -sf ${sampleFile}  \n\t -wd ${wkDir} \n\t -in ${initial} \n\t -set ${setFile} \n\t -db ${dbDir} \n"
${DIR}/manager.py -sf ${sampleFile} -wd ${wkDir} -in ${initial} -set ${setFile} -db ${dbDir}

echo ""
echo "---------------"
echo "MERGING RESULTS"
echo "---------------"
${DIR}/write_merged_xlsx.py -wd ${wkDir} -sf ${sampleFile} -in ${initial}

