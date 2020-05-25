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
  mkdir -p /data/famaprofiling/1.4
  cd /data/famaprofiling/1.4
  pwd

  echo "downloading nitrogen cycle database: http://iseq.lbl.gov/mydocs/fama_downloads/fama_nitrogen.tar.gz"
  curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/fama_nitrogen.tar.gz
  tar xvf /data/famaprofiling/1.4/fama_nitrogen.tar.gz
  rm /data/famaprofiling/1.4/fama_nitrogen.tar.gz
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.4/fama_nitrogen-cycle_classification_db_v.10.0.faa --db /data/famaprofiling/1.4/fama_nitrogen-cycle_classification_db_v.10.0
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.4/fama_nitrogen-cycle_preselection_db_v.10.0.faa --db /data/famaprofiling/1.4/fama_nitrogen-cycle_preselection_db_v.10.0
  rm /data/famaprofiling/1.4/fama_nitrogen-cycle_classification_db_v.10.0.faa
  rm /data/famaprofiling/1.4/fama_nitrogen-cycle_preselection_db_v.10.0.faa

  echo "downloading universal markers database: http://iseq.lbl.gov/mydocs/fama_downloads/fama_universal.tar.gz"
  curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/fama_universal.tar.gz
  tar xvf /data/famaprofiling/1.4/fama_universal.tar.gz
  rm /data/famaprofiling/1.4/fama_universal.tar.gz
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.4/fama_universal_classification_db_v.1.4.faa --db /data/famaprofiling/1.4/fama_universal-cycle_classification_db_v.1.4
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.4/fama_universal_preselection_db_v.1.4.faa --db /data/famaprofiling/1.4/fama_universal-cycle_preselection_db_v.1.4
  rm /data/famaprofiling/1.4/fama_universal_classification_db_v.1.4.faa
  rm /data/famaprofiling/1.4/fama_universal_preselection_db_v.1.4.faa

  echo "downloading nitrogen cycle database: http://iseq.lbl.gov/mydocs/fama_downloads/fama_rpl6.tar.gz"
  curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/fama_rpl6.tar.gz
  tar xvf /data/famaprofiling/1.4/fama_rpl6.tar.gz
  rm /data/famaprofiling/1.4/fama_rpl6.tar.gz
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.4/fama_rpl6_classification_db_v.1.2.faa --db /data/famaprofiling/1.4/fama_rpl6_classification_db_v.1.2
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.4/fama_rpl6_preselection_db_v.1.2.faa --db /data/famaprofiling/1.4/fama_rpl6_preselection_db_v.1.2
  rm /data/famaprofiling/1.4/fama_rpl6_classification_db_v.1.2.faa
  rm /data/famaprofiling/1.4/fama_rpl6_preselection_db_v.1.2.faa

  echo "downloading taxonomy database: http://iseq.lbl.gov/mydocs/fama_downloads/fama_taxonomy.tar.gz"
  curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/fama_taxonomy.tar.gz
  tar xvf /data/famaprofiling/1.4/fama_taxonomy.tar.gz
  rm /data/famaprofiling/1.4/fama_taxonomy.tar.gz

  echo "downloading Microbe Census data: http://iseq.lbl.gov/mydocs/fama_downloads/microbecensus_data.tar.gz"
  curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/microbecensus_data.tar.gz
  tar xvf /data/famaprofiling/1.4/microbecensus_data.tar.gz
  rm /data/famaprofiling/1.4/microbecensus_data.tar.gz
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.4/seqs.fa --db /data/famaprofiling/1.4/seqs
  rm /data/famaprofiling/1.4/seqs.fa

  if [ -s "/data/famaprofiling/1.4/seqs.dmnd" ] ; then
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
