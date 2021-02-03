SET DllPath=S:\OneDrive\OneDrive - Microsoft\Projects\2017 DNA Kinetics\DNA Structure Study\FastqProcessor\FastqProcessor\bin\Release\netcoreapp2.2\FastqProcessor.dll
SET SeqPath=S:\Source\Repos\BiologyGit\DnaModelling\Case Studies\DNA Kinetics\Structure Study\Data\StructureStudySequences.csv
SET DataPath=S:\Data\tmp\
start "" dotnet "%DllPath%" "%SeqPath%" "%DataPath%High_USPD16098953-4_H255VCCX2_L4_1.fq"
start "" dotnet "%DllPath%" "%SeqPath%" "%DataPath%hr0_USPD16098953-5_H255VCCX2_L4_1.fq"
start "" dotnet "%DllPath%" "%SeqPath%" "%DataPath%hr2_USPD16098953-6_H255VCCX2_L4_1.fq"
start "" dotnet "%DllPath%" "%SeqPath%" "%DataPath%hr4_USPD16098953-7_H255VCCX2_L4_1.fq"
start "" dotnet "%DllPath%" "%SeqPath%" "%DataPath%hr6_USPD16098953-12_H255VCCX2_L4_1.fq"
start "" dotnet "%DllPath%" "%SeqPath%" "%DataPath%Low_USPD16098953-2_H255VCCX2_L4_1.fq"