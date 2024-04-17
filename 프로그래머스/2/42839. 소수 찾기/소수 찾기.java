import java.util.*;
class Solution {
    int answer = 0;
    char[] chars;
    boolean[] visited;
    public int solution(String numbers) {
        chars = numbers.toCharArray();
        visited = new boolean[chars.length];
        HashSet<Integer> set = new HashSet<>();
        for(int i = 0; i < chars.length; i++){
            visited[i] = true;
            dfs(set, Integer.parseInt("" + chars[i]));
            visited[i] = false;
        }
        return answer;
    }
    public void dfs(HashSet<Integer> set, int value){
        if(!set.contains(value)){
            if(isPrimeNumber(value)) answer++;
            set.add(value);
        }
        for(int i = 0; i < visited.length; i++){
            if(visited[i]) continue;
            visited[i] = true;
            dfs(set, Integer.parseInt(value + "" + chars[i]));
            visited[i] = false;
        }
    }
    public boolean isPrimeNumber(int number){
        if(number <= 1) return false;
        for(int i = 2; i <= (int)Math.sqrt(number); i++){
            if(number % i == 0) return false;
        }
        return true;
    }
}