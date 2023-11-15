# Get FINEMAP binary
wget http://www.christianbenner.com/finemap_v1.4.2_x86_64.tgz
tar -xzf finemap_v1.4.2_x86_64.tgz
rm finemap_v1.4.2_x86_64.tgz
mv finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64 finemap
rm -rf finemap_v1.4.2_x86_64

# Get Hadoop AWS library and Java AWS SDK
wget https://repo1.maven.org/maven2/org/apache/hadoop/hadoop-aws/3.0.0/hadoop-aws-3.0.0.jar
wget https://repo1.maven.org/maven2/com/amazonaws/aws-java-sdk-bundle/1.11.199/aws-java-sdk-bundle-1.11.199.jar
