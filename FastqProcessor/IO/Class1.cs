using System;
using System.Runtime.InteropServices;
using System.IO;

namespace IO
{
    struct FastqEntry
    {
        string s;
    }

    public class Class1
    {
        


      [DllImport("msvcrt.dll", EntryPoint = "memcpy", CallingConvention = CallingConvention.Cdecl, SetLastError = false)]
        public static extern IntPtr memcpy(IntPtr dest, IntPtr src, UIntPtr count);

        unsafe static string Read(string path)
        {
            using (var file = File.OpenRead(path))
            {
                int size = sizeof(FastqEntry);
                const int BLOCK_SIZE = 512; // process at a time
                byte[] buffer = new byte[BLOCK_SIZE * size];

                UIntPtr bufferLen = new UIntPtr((uint)buffer.Length);
                fixed (byte* bufferPtr = buffer)
                {
                    Fill(file, buffer, 0, 4);
                    int len = ((int*)bufferPtr)[0];

                    FastqEntry[] result = new FastqEntry[len];
                    fixed (st* dataPtr = result)
                    {
                        st* rawPtr = dataPtr;
                        IntPtr source = new IntPtr(bufferPtr);
                        while (len >= BLOCK_SIZE)
                        {
                            Fill(file, buffer, 0, buffer.Length);
                            memcpy(new IntPtr(rawPtr), source, bufferLen);
                            len -= BLOCK_SIZE;
                            rawPtr += BLOCK_SIZE;
                        }
                        if (len > 0)
                        {
                            Fill(file, buffer, 0, len * size);
                            memcpy(new IntPtr(rawPtr), source, new UIntPtr((uint)(len * size)));
                        }
                    }
                    return result;
                }
            }


        }
        static void Fill(Stream source, byte[] buffer, int offset, int count)
        {
            int read;
            while (count > 0 && (read = source.Read(buffer, offset, count)) > 0)
            {
                offset += read;
                count -= read;
            }
            if (count > 0) throw new EndOfStreamException();
        }    
     
    }
}
