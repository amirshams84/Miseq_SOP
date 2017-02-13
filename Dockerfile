FROM amirshams/centos7:1.0

MAINTAINER Amir Shams <amir.shams84@gmail.com>
ENV ROOT=/
ENV CURRENT_PATH=.
CMD ["/bin/bash"]
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
RUN chmod -R 0755 $CURRENT_PATH/MISEQ_SOP_EXECDIR
RUN rm -rf $CURRENT_PATH/__MACOSX
RUN rm -rf $CURRENT_PATH/Mothur.linux_64.zip
##############################################################
# Software:             SILVA
# Software Version:     98
# Software Website:     -
# Description:          SILVA TAXONOMY FILE
##############################################################
RUN wget https://www.mothur.org/w/images/9/98/Silva.bacteria.zip
RUN unzip Silva.bacteria.zip
RUN mkdir $CURRENT_PATH/MISEQ_SOP_TAXDIR
RUN mv $CURRENT_PATH/silva.bacteria/* $CURRENT_PATH/MISEQ_SOP_TAXDIR
RUN chmod -R 0755 $CURRENT_PATH/MISEQ_SOP_TAXDIR
RUN rm -rf $CURRENT_PATH/__MACOSX
RUN rm -rf $CURRENT_PATH/silva.bacteria
RUN rm -rf $CURRENT_PATH/Silva.bacteria.zip
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
# Software:             PYTHON SCRIPT
# Software Version:     1
# Software Website:     -
# Description:          python
##############################################################
VOLUME $CURRENT_PATH/MISEQ_SOP_OUTPUTDIR
RUN wget https://raw.githubusercontent.com/amirshams84/Miseq_SOP/master/MiSeq_SOP.py -P $CURRENT_PATH/



CMD ["bin/bash"]