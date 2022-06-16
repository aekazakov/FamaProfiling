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
  mkdir -p /data/famaprofiling/1.5
  cd /data/famaprofiling/1.5
  pwd

  echo "downloading protein reference databases: http://iseq.lbl.gov/mydocs/fama_downloads/fama_nitrogen.tar.gz"
  curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/cgcms_fama_ref.tar.gz
  tar xvf /data/famaprofiling/1.5/cgcms_fama_ref.tar.gz
  rm /data/famaprofiling/1.5/cgcms_fama_ref.tar.gz
  cd /data/famaprofiling/1.5/fama/nitrogen11/
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.5/fama/nitrogen11/preselection_db_nr100.faa --db /data/famaprofiling/1.5/fama/nitrogen11/preselection_db_nr100
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.5/fama/nitrogen11/classification_db_nr100.faa --db /data/famaprofiling/1.5/fama/nitrogen11/classification_db_nr100
  rm /data/famaprofiling/1.5/fama/nitrogen11/preselection_db_nr100.faa
  rm /data/famaprofiling/1.5/fama/nitrogen11/classification_db_nr100.faa
  cd /data/famaprofiling/1.5/fama/universal1.4/
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.5/fama/universal1.4/preselection_db_nr100.faa --db /data/famaprofiling/1.5/fama/universal1.4/preselection_db_nr100
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.5/fama/universal1.4/classification_db_nr100.faa --db /data/famaprofiling/1.5/fama/universal1.4/classification_db_nr100
  rm /data/famaprofiling/1.5/fama/universal1.4/preselection_db_nr100.faa
  rm /data/famaprofiling/1.5/fama/universal1.4/classification_db_nr100.faa
  cd /data/famaprofiling/1.5/fama/cazy2/
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.5/fama/cazy2/preselection_database.faa --db /data/famaprofiling/1.5/fama/cazy2/preselection_database
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.5/fama/cazy2/classification_database.faa --db /data/famaprofiling/1.5/fama/cazy2/classification_database
  rm /data/famaprofiling/1.5/fama/cazy2/preselection_database.faa
  rm /data/famaprofiling/1.5/fama/cazy2/classification_database.faa
  
  echo "downloading ribosomal protein L6 database: http://iseq.lbl.gov/mydocs/fama_downloads/fama_rpl6.tar.gz"
  curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/fama_rpl6.tar.gz
  tar xvf /data/famaprofiling/1.5/fama_rpl6.tar.gz
  rm /data/famaprofiling/1.5/fama_rpl6.tar.gz
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.5/fama_rpl6_classification_db_v.1.2.faa --db /data/famaprofiling/1.5/fama_rpl6_classification_db_v.1.2
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.5/fama_rpl6_preselection_db_v.1.2.faa --db /data/famaprofiling/1.5/fama_rpl6_preselection_db_v.1.2
  rm /data/famaprofiling/1.5/fama_rpl6_classification_db_v.1.2.faa
  rm /data/famaprofiling/1.5/fama_rpl6_preselection_db_v.1.2.faa

  echo "downloading taxonomy database: http://iseq.lbl.gov/mydocs/fama_downloads/fama_taxonomy.tar.gz"
  curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/fama_taxonomy.tar.gz
  tar xvf /data/famaprofiling/1.5/fama_taxonomy.tar.gz
  rm /data/famaprofiling/1.5/fama_taxonomy.tar.gz

  echo "downloading Microbe Census data: http://iseq.lbl.gov/mydocs/fama_downloads/microbecensus_data.tar.gz"
  curl -LJO -q http://iseq.lbl.gov/mydocs/fama_downloads/microbecensus_data.tar.gz
  tar xvf /data/famaprofiling/1.5/microbecensus_data.tar.gz
  rm /data/famaprofiling/1.5/microbecensus_data.tar.gz
  /kb/deployment/bin/diamond/diamond makedb --in /data/famaprofiling/1.5/seqs.fa --db /data/famaprofiling/1.5/seqs
  rm /data/famaprofiling/1.5/seqs.fa

  if [ -s "/data/famaprofiling/1.5/seqs.dmnd" ] ; then
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
