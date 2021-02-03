module Dep
//open System.Runtime.InteropServices
//open System

//[<DllImport("msvcrt.dll", EntryPoint = "memcpy", CallingConvention = CallingConvention.Cdecl, SetLastError = false)>]
//extern IntPtr memcpy(IntPtr dest, IntPtr src, UIntPtr count);

//let Read(path:string) = 
//    use file = IO.File.OpenRead(path)
//    let size = sizeof<String>
//    let BLOCK_SIZE = 512; // process at a time
//    let buffer = Array.init (BLOCK_SIZE * size) (fun _ -> byte)

//    ""

////{
////    using (var file = File.OpenRead(path))
////    {//       

////        UIntPtr bufferLen = new UIntPtr((uint)buffer.Length);
////        fixed (byte* bufferPtr = buffer)
////        {
////            Fill(file, buffer, 0, 4);
////            int len = ((int*)bufferPtr)[0];

////            st[] result = new st[len];
////            fixed (st* dataPtr = result)
////            {
////                st* rawPtr = dataPtr;
////                IntPtr source= new IntPtr(bufferPtr);
////                while (len >= BLOCK_SIZE)
////                {
////                    Fill(file, buffer, 0, buffer.Length);
////                    memcpy(new IntPtr(rawPtr), source, bufferLen);
////                    len -= BLOCK_SIZE;
////                    rawPtr += BLOCK_SIZE;
////                }
////                if (len > 0)
////                {
////                    Fill(file, buffer, 0, len * size);
////                    memcpy(new IntPtr(rawPtr), source, new UIntPtr((uint)(len * size)));
////                }
////            }
////            return result;
////        }
////    }



//////let LoadFastq (file:string) = 
//////    file 
//////    |> System.IO.File.ReadAllLines
////    |> Array.chunkBySize 4
////    |> Array.map(fun block -> 
////        let sequence = block.[1]
////        let quality = block.[3]        
////        sequence, quality
////        )
    

////let rec GenerateFingerprints (k:int) (alphabet:string[]) = 
////    if k <= 0 then 
////        alphabet |> Seq.toArray |> Array.map string 
////    else
////        GenerateFingerprints (k-1) alphabet                
////        |> Array.collect(fun (s:string) -> 
////            alphabet
////            |> Seq.toArray
////            |> Array.map(fun c -> s + (string c))
////            )        


////let labelFingerprint (s:string) (F:string[]) = F |> Array.indexed |> Array.filter(fun (_,f) -> s.Contains f ) |> Array.map fst |> Set.ofSeq
    
    




////let LabelReadsBasic (sequences:(string*string)[]) (reads:(string*string)[]) = 
////        let t0 = System.DateTime.Now        
////        let S = sequences |> Array.map (fun (name,s) -> name, s.Replace("T","*").Replace("C","*"), s)
        
////        let labelled = 
////            reads              
////            |> Array.choose(fun (s0,_) -> 
////                let s = s0.Replace("T","*").Replace("C","*")                                        
////                //let full =
////                let full = S |> Array.tryFind (fun (_,s',_) -> s.Contains s')
////                match full with 
////                | Some (name,s',s'') ->                     
////                    //let i = s.IndexOf s'
////                    //printfn "\n%s\n%s\n%s\n%s%s" (if s0.Contains s'' then "exact" else "approx") name s0 (String.replicate i ".") s''
////                    Some name
////                | None -> None        
////            )            
////        printfn "%i reads total" reads.Length
////        printfn "%i reads labelled" labelled.Length
////        printfn "%f of the data was used" ((float labelled.Length) / (float reads.Length))
////        printfn "%.2f seconds"  (System.DateTime.Now - t0).TotalSeconds

////let LabelReadsAho (sequences:(string*string)[]) (reads:(string*string)[]) = 
////        let t0 = System.DateTime.Now        
////        let S = 
////            sequences 
////            |> Array.map (fun (name,s) -> name, s.Replace("T","*").Replace("C","*"), s)
////            |> Array.map(fun (_,s,_) -> s)
////            |> AhoCorasick
        
////        let labelled = 
////            reads              
////            |> Array.Parallel.choose(fun (s0,_) -> 
////                let s = s0.Replace("T","*").Replace("C","*")                                        
////                if s |> S.Search |> Seq.isEmpty then None else Some ""
////            )            
////        printfn "%i reads total" reads.Length
////        printfn "%i reads labelled" labelled.Length
////        printfn "%f of the data was used" ((float labelled.Length) / (float reads.Length))
////        printfn "%.2f seconds"  (System.DateTime.Now - t0).TotalSeconds
     

     

////let LabelReadsMerged (sequences:(string*string)[]) (reads:(string*string)[]) = 
////        let t0 = System.DateTime.Now        
////        let Smap = 
////            sequences 
////            |> Array.map (fun (name,s) -> 
////                let s' = s.Replace("T","*").Replace("C","*")
////                name, s', s.Length)
////        let S = Smap |> Array.map(fun (_,s,_) -> s) |> String.concat ""

////        let labelled = 
////            reads              
////            |> Array.choose(fun (s0,_) -> 
////                let s = s0.Replace("T","*").Replace("C","*")                                                        
////                if S.IndexOf(s)>0 then                 
////                    Some ""                    
////                else None
////                )            
////        printfn "%i reads total" reads.Length
////        printfn "%i reads labelled" labelled.Length
////        printfn "%f of the data was used" ((float labelled.Length) / (float reads.Length))
////        printfn "%.2f seconds"  (System.DateTime.Now - t0).TotalSeconds
     
////let LabelReadsFingerprint k (sequences:(string*string)[]) (reads:(string*string)[]) = 
////        let t0 = System.DateTime.Now        
////        let F = GenerateFingerprints k [|"A"; "G"; "*"|]               
////        let S = 
////            sequences 
////            |> Array.map (fun (name,s) -> 
////                let s' = s.Replace("T","*").Replace("C","*")
////                let label = labelFingerprint s' F
////                label, (name, s', s))            
        
////        let labelled = 
////            reads              
////            |> Array.choose(fun (s0,_) -> 
////                let s = s0.Replace("T","*").Replace("C","*")                                        
////                let l = labelFingerprint s F
                                    
////                let S' = S |> Array.filter(fun (L,_) -> Set.isSubset L l) |> Array.map snd
                
////                let full =  S' |> Array.tryFind (fun (_,s',_) -> s.Contains s')
////                match full with 
////                | Some (name,s',s'') ->                     
////                    //let i = s.IndexOf s'
////                    //printfn "\n%s\n%s\n%s\n%s%s" (if s0.Contains s'' then "exact" else "approx") name s0 (String.replicate i ".") s''
////                    Some name
////                | None -> None  
                
////            )            
////        printfn "%i reads total" reads.Length
////        printfn "%i reads labelled" labelled.Length
////        printfn "%f of the data was used" ((float labelled.Length) / (float reads.Length))
////        printfn "%.2f seconds"  (System.DateTime.Now - t0).TotalSeconds

////let LabelReadsAlign (sequences:(string*string)[]) (reads:(string*string)[]) = 
////        let t0 = System.DateTime.Now        
////        let toSeq (s:string) = Bio.Sequence(Bio.Alphabets.DNA, s)

////        let S = sequences |> Array.map(fun (name,s) -> name, toSeq s)

////        //let aligner = Bio.Algorithms.Alignment.NeedlemanWunschAligner()            
////        let aligner = Bio.Algorithms.Alignment.SmithWatermanAligner()        
////        //aligner.GapOpenCost <- -1
////        //aligner.GapExtensionCost <- 0
////        //let labelled =                         
////        reads          
////        |> Array.map(fun (s0,_) ->             
////            let s = toSeq s0                
////            //let flag = sequences |> Array.tryFind(fun (name,s') -> (s0.Replace("T","*").Replace("C","*")).Contains (s'.Replace("T","*").Replace("C","*")))
////            //printfn "%A" flag                        
////            S
////            |> Array.map(fun (name,s') ->                     
////                let al = aligner.AlignSimple(s,s').[0]                
////                name, s', al.PairwiseAlignedSequences.[0]
////                )
////            |> Array.sortByDescending (fun (_,_,x) -> x.Score)//x.FirstSequence.Count)
////            |> fun X -> X.[0]
////            )
////        |> Array.iteri(fun i (name, s', x) ->            
////            printfn "Read %i: %i, %i" i (s'.Count - x.FirstSequence.Count) x.Score
////            )

////            //|> Array.take 5
////            //|> Array.iter(fun (name,score,al) -> 
////            //    //al.ToString() |> printfn "%s"

////            //    //let p1 = String.replicate (int al.FirstOffset) "."
////            //    //let p2 = String.replicate (int al.SecondOffset) "."                
////            //    //printfn "%s (%i)\n%s%s\n%s%s" name score p2 (al.FirstSequence.ToString()) p1 (al.SecondSequence.ToString())
////            //    )
                                   
            
////        printfn "%i reads total" reads.Length
////        //printfn "%i reads labelled" labelled.Length
////        //printfn "%f of the data was used" ((float labelled.Length) / (float reads.Length))
////        printfn "%.2f seconds"  (System.DateTime.Now - t0).TotalSeconds
        

  
     
////let LabelReadsAlignFingerprint k (sequences:(string*string)[]) (reads:(string*string)[]) = 
////        let t0 = System.DateTime.Now        
////        let toSeq (s:string) = Bio.Sequence(Bio.Alphabets.DNA, s)
////        let F = GenerateFingerprints k [|"A"; "C"; "G"; "T"|]        
////        let S = 
////            sequences 
////            |> Array.map(fun (name,s) -> name, labelFingerprint s F, toSeq s)                        
            
////        //let aligner = Bio.Algorithms.Alignment.NeedlemanWunschAligner()            
////        let aligner = Bio.Algorithms.Alignment.SmithWatermanAligner()        
////        //aligner.GapOpenCost <- -1
////        //aligner.GapExtensionCost <- 0
////        //let labelled =                         
////        reads          
////        |> Array.indexed
////        |> Array.choose(fun (i,(s0,_)) ->             
////            let s = toSeq s0                
////            //let flag = sequences |> Array.tryFind(fun (name,s') -> (s0.Replace("T","*").Replace("C","*")).Contains (s'.Replace("T","*").Replace("C","*")))
////            //printfn "%A" flag                        
////            let f = labelFingerprint s0 F
////            let nt = (int (0.4*(float f.Count)))
            
////            let S' = S |> Array.filter(fun (_,f',_) -> (Set.intersect f f').Count > nt)

////            if Array.isEmpty S' then 
////                None
////            else
////                S'
////                |> Array.map(fun (name,_,s') ->                 
////                    let al = aligner.AlignSimple(s,s')                             
////                    name, s', al.[0].PairwiseAlignedSequences.[0]
////                    )
////                |> Array.sortByDescending (fun (_,_,x) -> x.Score)//x.FirstSequence.Count)
////                |> fun X -> Some (i, X.[0])
            
////            )
////        |> Array.iter(fun (i, (name, s', x)) ->            
////            printfn "Read %i: %i, %i" i (s'.Count - x.FirstSequence.Count) x.Score
////            x.ToString() |> printfn "%s\n"
////            )

////            //|> Array.take 5
////            //|> Array.iter(fun (name,score,al) -> 
////            //    //al.ToString() |> printfn "%s"

////            //    //let p1 = String.replicate (int al.FirstOffset) "."
////            //    //let p2 = String.replicate (int al.SecondOffset) "."                
////            //    //printfn "%s (%i)\n%s%s\n%s%s" name score p2 (al.FirstSequence.ToString()) p1 (al.SecondSequence.ToString())
////            //    )
                                   
            
////        printfn "%i reads total" reads.Length
////        //printfn "%i reads labelled" labelled.Length
////        //printfn "%f of the data was used" ((float labelled.Length) / (float reads.Length))
////        printfn "%.2f seconds"  (System.DateTime.Now - t0).TotalSeconds
