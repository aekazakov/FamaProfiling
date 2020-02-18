FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update && \
    apt-get -y install automake build-essential zlib1g libz-dev

# Install DIAMOND 

RUN mkdir -p /kb/deployment/bin/diamond && \
    git clone https://github.com/bbuchfink/diamond.git && \
    cd diamond && \
    git checkout e186f17 . && \
    sed 's/g++/g++ -std=c++11/' build_simple.sh > tmp && \
    chmod u+x tmp && \
    mv tmp rebuild_simple.sh && \
    ./rebuild_simple.sh && \
    mv diamond /kb/deployment/bin/diamond

RUN /kb/deployment/bin/diamond/diamond version

# Install Krona Tools
RUN mkdir -p /kb/deployment/bin/krona && \
    cd /kb/deployment/bin/krona && \
    git clone https://github.com/marbl/Krona && \
    cd Krona/KronaTools && \
    ./install.pl --prefix /kb/deployment/bin/

RUN pip install statsmodels mathstats xlsxwriter fpdf

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

# Prepare directories for DIAMOND databases (rest of db installation to ref data mount in entrypoint.sh init script)
RUN mkdir -p /data/famaprofiling

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
