class Solution {
    public String solution(String polynomial) {
        String[] sp = polynomial.split(" [+] ");
        int x = 0;
        int n = 0;
        for(int i = 0; i < sp.length; i++){
            if(sp[i].equals("x")) x++;
            else if(sp[i].contains("x")) x += Integer.parseInt(sp[i].replaceAll("x", ""));
            else n += Integer.parseInt(sp[i]);
        }
        String answer = "";
        if(x > 1) answer += x;
        if(x > 0) answer += "x";
        if(x > 0 && n > 0) answer += " + ";
        if(n > 0) answer += n;
        return answer;
        
    }
}