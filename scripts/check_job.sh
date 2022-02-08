#!/bin/bash

set -x

# script to print message and cancel if jobs never satisfied or cancelled
# change : to ;
jids=${jid_to_check//:/,}
# states of jobs in queue
state_q=($(squeue -j $jids -h))

# while jobs in queue keep monitoring
while [ ${#state_q[@]} -ne 0 ];do
        # if any job  has dependency that never satisfied or if they are cancelled or failed
        # cancel all queued jobs.
        reason=($(squeue -j $jids -o "%R" -h))
        state=($(sacct -j $jids --format=state | tail -n +3 ))
	# store jobname and jobids
	jobname=($(sacct -j $jids --format=jobname%100 | tail -n +3 ))
	jobid=($(sacct -j $jids --format=jobid | tail -n +3 ))
	# unique states of jids
	state_uni=($(echo "${state[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
        if [[ "${reason[*]}" == *"DependencyNeverSatisfied"* || "${state[*]}" == *"CANCELLED"* || \
                "${state[*]}" == *"FAILED"* || "${reason[*]}" == *"FAILED"* ]]; then
                for i in "${!state[@]}"; do
   			if [[ "${state[$i]}" = *"FAILED"* ]]; then
       				fail_idx=$i;
				break
   			fi
		done
		scancel "$jids"
                echo -e "$msg_fail" >> "${out_file}"
		fail_job_info="${jobname[$fail_idx]} jobid:${jobid[$fail_idx]} failed.\n"
		fail_job_info="${fail_job_info}See ${jobname[$fail_idx]}.out for more details.\n\n"
		echo -e "$fail_job_info" >> "${out_file}"
		if [[ -f "$(dirname $out_file)/run_rnaseq_full.out" ]];then
			echo -e "$fail_job_info" >> "$(dirname $out_file)/run_rnaseq_full.out"
		fi
                cp "${out_file}" "$(dirname $out_file)/outputs/logs/$(basename $out_file)"
                date >> "${out_file}"
                exit 1
	elif [[ "${#state_uni[@]}" -eq 1 ]] && [[ "${state_uni[0]}" == *"COMPLETED"* ]] ;then
                # jobs completed exit early
		break
        else
                sleep 10
        fi
	# update states of jobs in queue
	state_q=($(squeue -j $jids -h))
	echo $state_q
done

# jobs should be successfully completed 
# print ok message
echo -e "$msg_ok" >> "${out_file}"
date >> "${out_file}"
echo -e "\n" >> "${out_file}"
cp "${out_file}" "$(dirname $out_file)/outputs/logs/$(basename $out_file)"
