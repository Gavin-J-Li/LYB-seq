// Learn more about F# at http://fsharp.org

open System
open Ganss.Text

let ProcessFastq f (path:string) =     
    let fs = System.IO.File.Open(path, System.IO.FileMode.Open, System.IO.FileAccess.Read, System.IO.FileShare.ReadWrite)
    let bs = new System.IO.BufferedStream(fs)
    let sr = new System.IO.StreamReader(bs)    
    
    while not sr.EndOfStream do        
        let seqId = sr.ReadLine()
        let s = sr.ReadLine()
        sr.ReadLine() |> ignore
        sr.ReadLine() |> ignore        
        f seqId s    

        

type Conversion = 
    { baseId : int
      count  : int
    }
    static member Create (i, x)  = {baseId = i; count = x}

type ConvResult = 
    { name : string
      correct_reads: int
      error_reads: int
      conv : Conversion[]
      error: Conversion[]
    }


[<EntryPoint>]
let main argv =    
    
    //let DataFolder = @"S:\Source\Repos\BiologyGit\DnaModelling\Case Studies\DNA Kinetics\Structure Study\Data"
    let sequence_file = argv.[0]
    let data_file = argv.[1]
    let dataset_name = data_file.[data_file.LastIndexOf("\\")+1..data_file.LastIndexOf(".")-1]

    printfn "Converting %s\n\n" dataset_name

    let sequences = 
        sequence_file
        |> System.IO.File.ReadAllLines
        |> Array.map(fun l -> 
            let f = l.Split(',')
            f.[0], f.[1])
        |> Array.map (fun (name,s) -> s.Replace("T","*").Replace("C","*"), s, name)
           
    let sequenceDNA = sequences |> Array.map(fun (_,s,_) -> s)

    let sequenceIndex = 
        sequences
        |> Array.map (fun (s,_,_) -> s)
        |> AhoCorasick
    
    let sequenceMap = 
        sequences 
        |> Array.mapi (fun i (s,s',name) -> s, (s',name, i))
        |> Map.ofSeq

    let t0 = System.DateTime.Now
        
    
    let mutable Ntotal = 0    
    let mutable Nerror = 0 //T->C conversion
    let mutable Ncorrect = 0        
    let mutable reads = List.empty

    let conv_rate  = Array.init sequences.Length (fun i -> Array.init sequenceDNA.[i].Length (fun _ -> 0))
    let error_rate = Array.init sequences.Length (fun i -> Array.init sequenceDNA.[i].Length (fun _ -> 0))
    let correct_reads = Array.init sequences.Length (fun _ -> 0)
    let error_reads = Array.init sequences.Length (fun _ -> 0)


    IO.Directory.CreateDirectory(dataset_name) |> ignore
    
    data_file
    |> ProcessFastq (fun readId s0 ->         
        let s = s0.Replace("T","*").Replace("C","*")        
        Ntotal <- Ntotal + 1
                    
        if Ntotal%100000=0 then 
            let num2str (x:int) = x |> sprintf "%i" |> Seq.rev |> Seq.chunkBySize 3 |> Seq.map System.String |> String.concat "," |> Seq.rev |> Array.ofSeq |> System.String
            printfn "%s reads processed\n\t%s errors\n\t%s correct\n\t%.2f of data was labelled\n" (num2str Ntotal) (num2str Nerror) (num2str Ncorrect) ((float Ncorrect)/(float Ntotal))
            printfn "%f lines per second" ((float Ntotal)/(System.DateTime.Now - t0).TotalSeconds)

        let labels = sequenceIndex.Search s |> Array.ofSeq
            
        if  labels.Length > 0 then                                           
            //readID: NGS read
            //name: name of template sequence
            //s': NGS read sequence (aligned)
            //s: template sequence (aligned)

            let m = labels.[0] //can there be multiple labels?                
            let s, name, seqId = sequenceMap.[m.Word]
            let j = s.Length + m.Index - 1
            let s' = s0.[m.Index..j]                                
            let S = Seq.zip s s'             

            //first pass, check for errors
            let mutable has_errors = false
            S
            |> Seq.iteri(fun i (c,c') -> 
                if c='T' && c'='C' then 
                    has_errors <- true
                    error_rate.[seqId].[i] <- error_rate.[seqId].[i] + 1                    
                )
            
            if has_errors then
                Nerror <- Nerror + 1          
                error_reads.[seqId] <- error_reads.[seqId] + 1
            else
                Ncorrect <- Ncorrect + 1                                        
                correct_reads.[seqId] <- correct_reads.[seqId] + 1                                    
                IO.File.AppendAllLines(sprintf "%s/%s.txt" dataset_name name, [|s'|])
                //reads <- (name,s')::reads
                S
                |> Seq.iteri(fun i (c,c') -> if c='C' && c'='T' then conv_rate.[seqId].[i] <- conv_rate.[seqId].[i] + 1)
                )     

    let conv_rate = conv_rate |> Array.map (Array.indexed >> Array.filter(fun (_,x) -> x > 0))
    let error_rate = error_rate |> Array.map (Array.indexed >> Array.filter(fun (_,x) -> x > 0))


    sequences
    |> Array.mapi(fun i (_,_,name) ->         
        { name  = name
          correct_reads = correct_reads.[i]
          error_reads = error_reads.[i]
          conv  = conv_rate.[i] |> Array.map Conversion.Create
          error = error_rate.[i] |> Array.map Conversion.Create
         }
        )
    |> Newtonsoft.Json.JsonConvert.SerializeObject
    |> fun x -> IO.File.WriteAllText(sprintf "%s/counts.json" dataset_name, x)
    
    [ sprintf "%i reads total" Ntotal
      sprintf "%i reads have errors" Nerror
      sprintf "%i reads labelled correctly" Ncorrect
      sprintf "%f of the data was used" ((float Ncorrect) / (float Ntotal))
      sprintf "%.2f seconds"  (System.DateTime.Now - t0).TotalSeconds
    ]
    |> String.concat "\n"
    |> fun x -> IO.File.WriteAllText(sprintf "%s/log.txt" dataset_name, x)
          
   
    //reads
    //|> Seq.groupBy fst
    //|> Seq.iter(fun (name, L) -> 
    //    IO.File.WriteAllLines(sprintf "%s/%s.txt" dataset_name name, L |> Seq.map snd)
    //    )
    0 // return an integer exit code
