FROM amirshams/centos7:1.0

MAINTAINER Amir Shams <amir.shams84@gmail.com>
ENV ROOT=/
ENV CURRENT_PATH=.
CMD ["/bin/bash"]
##############################################################
# Software:             PIP INSTALL PACKAGES
# Software Version:     1.0
# Software Website:     -
# Description:          required javascript library
##############################################################
RUN pip install numpy
RUN pip install scipy
RUN pip install plotly
RUN pip install pandas
RUN pip install biom-format
RUN pip install xlrd
RUN pip install openpyxl
RUN pip install xlwt
RUN pip install XlsxWriter
RUN pip install lxml
RUN pip install zip
##############################################################
# Software:             MOTHUR
# Software Version:     98
# Software Website:     -
# Description:          SILVA TAXONOMY FILE
##############################################################
RUN wget https://github.com/mothur/mothur/releases/download/v1.39.1/Mothur.linux_64.zip
RUN unzip Mothur.linux_64.zip
RUN mkdir $CURRENT_PATH/MISEQ_SOP_EXECDIR
RUN mv $CURRENT_PATH/mothur $CURRENT_PATH/MISEQ_SOP_EXECDIR/
RUN rm -rf $CURRENT_PATH/__MACOSX
RUN rm -rf $CURRENT_PATH/Mothur.linux_64.zip
##############################################################
# Software:             MAFFT
# Software Version:     7.22
# Software Website:     -
# Description:          MAFFT
##############################################################
RUN wget http://mafft.cbrc.jp/alignment/software/mafft-7.222-linux.tgz
RUN tar zxvf mafft-7.222-linux.tgz
RUN mv $CURRENT_PATH/mafft-linux64/* $CURRENT_PATH/MISEQ_SOP_EXECDIR/
RUN rm -rf $CURRENT_PATH/mafft-7.222-linux.tgz
RUN rm -rf $CURRENT_PATH/mafft-linux64
RUN rm -rf $CURRENT_PATH/mafft-linux32
RUN chmod -R 0755 $CURRENT_PATH/MISEQ_SOP_EXECDIR
##############################################################
# Software:             TEST DATA
# Software Version:     1
# Software Website:     -
# Description:          TEST DATA
##############################################################
RUN wget https://www.mothur.org/w/images/d/d6/MiSeqSOPData.zip
RUN unzip MiSeqSOPData.zip
RUN mkdir $CURRENT_PATH/MISEQ_SOP_TESTDIR
RUN mv $CURRENT_PATH/MiSeq_SOP/* $CURRENT_PATH/MISEQ_SOP_TESTDIR
RUN chmod -R 0755 $CURRENT_PATH/MISEQ_SOP_TESTDIR
RUN rm -rf $CURRENT_PATH/__MACOSX
RUN rm -rf $CURRENT_PATH/MiSeq_SOP
RUN rm -rf $CURRENT_PATH/MiSeqSOPData.zip
##############################################################
# Software:             GREENGENESE
# Software Version:     13_8_99
# Software Website:     -
# Description:          gg_13_8_99
##############################################################
RUN wget http://www.mothur.org/w/images/6/68/Gg_13_8_99.taxonomy.tgz
RUN tar zxvf Gg_13_8_99.taxonomy.tgz
RUN mkdir $CURRENT_PATH/MISEQ_SOP_TAXDIR
RUN mv $CURRENT_PATH/gg_13_8_99* $CURRENT_PATH/MISEQ_SOP_TAXDIR
RUN rm -rf $CURRENT_PATH/README.Rmd
RUN rm -rf $CURRENT_PATH/README.html
RUN rm -rf $CURRENT_PATH/Gg_13_8_99.taxonomy.tgz
RUN chmod -R 0755 $CURRENT_PATH/MISEQ_SOP_TAXDIR
##############################################################
# Software:             PYTHON SCRIPT
# Software Version:     1
# Software Website:     -
# Description:          python
##############################################################
VOLUME $CURRENT_PATH/MISEQ_SOP_OUTPUTDIR
VOLUME $CURRENT_PATH/MISEQ_SOP_INPUTDIR
RUN wget https://raw.githubusercontent.com/amirshams84/Miseq_SOP/master/miseq_sop.py -P $CURRENT_PATH/



CMD ["bin/bash"]