using System;
using System.Collections;
public class Solution {
    public int solution(string[,] clothes) {
        Hashtable hash = new Hashtable();
        for(int i = 0; i < clothes.GetLength(0); i++){
            string key = clothes[i,1];
            if(hash.ContainsKey(key)){
                int val = (int)hash[key];
                hash.Remove(key);
                hash.Add(key, val + 1);
            }
            else hash.Add(key, 1);
        }
        int count = 1;
        foreach(int v in hash.Values){
            count *= v + 1;
        }
        return count - 1;
    }
    int factorial(int n){
        if(n == 1) return 1;
        return n * factorial(n - 1);
    }
}