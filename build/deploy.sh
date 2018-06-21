

#This variable must be iqual to id in distributionManagement.repository/snapshotRepository on pom.xml
DEPLOYMENT_REPO_ID="oss.sonatype.org"

#echo OSS_SONATYPE_USER $OSS_SONATYPE_USER
#echo OSS_SONATYPE_PASS $OSS_SONATYPE_PASS



H2_HOME=$HOME/.deploym2

echo ${H2_HOME}

M2REPOSITORY=$H2_HOME/repository
mkdir -p ${M2REPOSITORY}
cp build/oss.sonatype.settings.xml ${H2_HOME}/settings.xml



ls ${H2_HOME}

sed -i -e s,LOCAL_M2,$H2_HOME,g ${H2_HOME}/settings.xml
sed -i -e s/OSS_SONATYPE_USER/$OSS_SONATYPE_USER/g ${H2_HOME}/settings.xml
sed -i -e s/OSS_SONATYPE_PASS/$OSS_SONATYPE_PASS/g ${H2_HOME}/settings.xml
sed -i -e s/DEPLOYMENT_REPO_ID/$DEPLOYMENT_REPO_ID/g ${H2_HOME}/settings.xml

#GPG_PRIVATE_KEY=`cat build/gpg_private`
#echo GPG_PRIVATE_KEY: $GPG_PRIVATE_KEY
#gpg -v --batch --import <(echo "$GPG_PRIVATE_KEY")

mkdir -p $HOME/lib/
sshpass -p ${SSH_PASS} scp -v -o "StrictHostKeyChecking no" ${SSH_USER}@${SERVER_WITH_CPLEX}:cplex.jar $HOME/lib/cplex.jar
ls -lh $HOME/
ls -lh $HOME/lib

mvn -DskipTests=true -Dcplex.jar.path=$HOME/lib/cplex.jar -s ${H2_HOME}/settings.xml -Dmaven.repo.local=$M2REPOSITORY deploy