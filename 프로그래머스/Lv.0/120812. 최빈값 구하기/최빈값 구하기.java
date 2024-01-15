import java.util.*;
class Max{
    int count;
    int value;
    public Max(int count, int value){
        this.count = count;
        this.value = value;
    }
}
class Solution {
    public int solution(int[] array) {
        Arrays.sort(array);
        int count = 0;
        int value = -1;
        List<Max> list = new ArrayList<>();
        for(int i = 0; i < array.length; i++){
            if(value != array[i]){
                if(count > 0){
                    if(list.size() == 0 || list.get(0).count <= count){
                        if(list.size() != 0 && list.get(0).count < count) list.clear();
                        list.add(new Max(count, value));
                    }
                }
                count = 1;
                value = array[i];
            }else count++;
            
            if(i == array.length - 1){
                if(list.size() == 0 || list.get(0).count <= count){
                    if(list.size() != 0 && list.get(0).count < count) list.clear();
                    list.add(new Max(count, value));
                }
            }
        }
        return list.size() > 1 ? -1 : list.get(0).value;
    }
}