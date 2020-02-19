#!/bin/bash

. /kb/deployment/user-env.sh

python ./scripts/prepare_deploy_cfg.py ./deploy.cfg ./work/config.properties

if [ -f ./work/token ] ; then
  export KB_AUTH_TOKEN=$(<./work/token)
fi

if [ $# -eq 0 ] ; then
  sh ./scripts/start_server.sh
elif [ "${1}" = "test" ] ; then
  echo "Run Tests"
  make test
elif [ "${1}" = "async" ] ; then
  sh ./scripts/run_async.sh
elif [ "${1}" = "init" ] ; then
  echo "Initialize module"
  mkdir -p /data/famaprofiling/1.3
  cd /data/famaprofiling/1.3
  pwd

  echo "downloading nitrogen cycle database: http://iseq.lbl.gov/mydocs/fama_downloads/fama_nitrogen_v10.tar.gz"
  curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/fama_nitrogen_v10.tar.gz
  tar xvf /data/famaprofiling/1.3/fama_nitrogen_v10.tar.gz
  rm /data/famaprofiling/1.3/fama_nitrogen_v10.tar.gz
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.3/fama_nitrogen-cycle_classification_db_v.10.0.faa --db /data/famaprofiling/1.3/fama_nitrogen-cycle_classification_db_v.10.0
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.3/fama_nitrogen-cycle_preselection_db_v.10.0.faa --db /data/famaprofiling/1.3/fama_nitrogen-cycle_preselection_db_v.10.0
  rm /data/famaprofiling/1.3/fama_nitrogen-cycle_classification_db_v.10.0.faa
  rm /data/famaprofiling/1.3/fama_nitrogen-cycle_preselection_db_v.10.0.faa

  echo "downloading taxonomy database: http://iseq.lbl.gov/mydocs/fama_downloads/fama1_taxonomy.tar.gz"
  curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/fama1_taxonomy.tar.gz
  tar xvf /data/famaprofiling/1.3/fama1_taxonomy.tar.gz
  rm /data/famaprofiling/1.3/fama1_taxonomy.tar.gz

  echo "downloading Microbe Census data: http://iseq.lbl.gov/mydocs/fama_downloads/microbecensus_data.tar.gz"
  curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/microbecensus_data.tar.gz
  tar xvf /data/famaprofiling/1.3/microbecensus_data.tar.gz
  rm /data/famaprofiling/1.3/microbecensus_data.tar.gz
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.3/seqs.fa --db /data/famaprofiling/1.3/seqs
  rm /data/famaprofiling/1.3/seqs.fa

  if [ -s "/data/famaprofiling/seqs.dmnd" ] ; then
    echo "DATA DOWNLOADED SUCCESSFULLY"
    touch /data/__READY__
  else
    echo "Init failed"
  fi
elif [ "${1}" = "bash" ] ; then
  bash
elif [ "${1}" = "report" ] ; then
  export KB_SDK_COMPILE_REPORT_FILE=./work/compile_report.json
  make compile
else
  echo Unknown
fi
