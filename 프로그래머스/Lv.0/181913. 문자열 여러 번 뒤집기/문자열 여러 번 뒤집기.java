class Solution {
    public String solution(String my_string, int[][] queries) {
        char[] chars = my_string.toCharArray();
        for(int i = 0; i < queries.length; i++){
            while(queries[i][0] < queries[i][1]){
                char tmp = chars[queries[i][0]];
                chars[queries[i][0]] = chars[queries[i][1]];
                chars[queries[i][1]] = tmp;
                queries[i][0]++;
                queries[i][1]--;
            }
        }
        return String.valueOf(chars);
    }
}