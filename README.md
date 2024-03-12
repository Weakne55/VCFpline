### Для работы с пайплайном необходимо сделать несколько действий перед началом работы:
1) Создать рабочую директорию и скачать в неё репозиторий
```
mkdir pipline
cd pipline
git clone https://github.com/Weakne55/VCFpline.git
```
2) Скачать micromamba и установить зависимости:
```
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
micromamba create -f env.yml
```
После этого необходимо активировать это окружение
```
micromamba activate vcf
```
Теперь слева или справа от названия оболочки вашей командной строки должна появится надпись **vcf** в скобках
3) Теперь нужно скачать геном человека, для последующей работы с ним
```
datasets download genome accession GCF_000001405.26
```
Распакуем его 
```
unzip ncbi_dataset.zip
rm ncbi_dataset.zip
```
### Начало работы
Для работы пайплайна нужно перенести в директорию **pipline** папку, в которой содержатся файлы секвенирования вида name_R[1,2]_001.fastq.gz 
Чтобы запустить пайплайн напишите команду
```
bash SNPline <folder_name>
```
или
```
./SNPline <folder_name>
```
### Результаты
В результате работы будет создана папка snpline_results внутри которой будут содержаться поддериктории согласно названиям файлов содержащихся в исходной папке с данными.
Внтури каждой такой директории будут находится 
Папки: 
    - FASTQC_POST,FASTQC_PRE - в которых содержатся результаты работы fastqc до и после удаления адаптеров с ридов
    - adapterremoval - в которой содержатся данные после удаления адаптеров, т.е. работы программы AdapterRemoval
    - results_aln_mapped_sorted_rmdup -  в которой содержаться результаты работы программы MapDamage
Файлы:
    - Исходные данные в формате name_R[1,2]_001.fastq.gz
    - aln_mapped_sorted.bam - риды, которые были выровнены на геном с помощью bwa aln и отсортированы с помощью samtools
    - aln_mapped_sorted_rmdup.bam - риды без дупликатов
    - aln_metrics.txt - статистика работы MarkDuplicates
    - bcfmpileup.bcf - файл с правдободобностью нахождения замен в геноме (liklihoods)
    - bcfview.bcf - отсортированные данные из bcfmpileup.bcf для поиска SNP
    - bcfcall.vcf - итговый файл, содержащий SNP 
### Примеры
SNP на митохондриальном геноме
![мито](https://github.com/Weakne55/VCFpline/assets/72615087/b684ddaf-2c95-4772-84f8-698a1233e7fe) 
SNP на Y хромосоме 
![игрек](https://github.com/Weakne55/VCFpline/assets/72615087/1d70723b-cc39-4f5a-80ed-35f1845ad2ba)
