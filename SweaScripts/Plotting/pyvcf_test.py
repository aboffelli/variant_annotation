import vcf

with open('/Users/student/Box/Notes/TestData/CustomAnnotation/custom_edit_gene_'
          'name_binding_sites_edited_fs_filtered_Lu0002_L05319_A19_S3.r1_VEP.'
          'vcf', 'r') as file:
    vcf_file = file.readlines()

vcf_test = vcf.Reader(vcf_file)

for line in vcf_test:
    print(line.INFO['CSQ'])
