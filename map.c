    /*
    Written by William Sutherland for 
    COMP20007 Assignment 1 2023 Semester 1
    Modified by Grady Fitzpatrick
    
    Implementation for module which contains map-related 
        data structures and functions.

    Functions in the task description to implement can
        be found here.
    
    Code implemented by <YOU>
*/
#include "map.h"
#include "stack.h"
#include "pq.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
int test_coord(int x_coord, int y_coord, int height, int width);
int comparator(const void *p, const void *q);
struct point ** getgraph(struct map *m);
void dfs(int* visited, struct point** adj_list, int no, int* num_nodes, int* component_sum, struct map *m);
int get_value(int x, int y, struct map *m);
int get_node_number(int x, int y, struct map *m);
struct point** add_airports(struct point ** adj_list, struct map *m);
int djikstra(struct map *m, struct point *start, struct point *end, struct point ** adj_list);
double get_time(struct point *start, struct point *end, struct point** adj_list, struct map *m);
struct point* get_closest_airport(struct point *start, struct point** adj_list, int num_points, struct map *m, int* time);
int get_num_gaps(struct point *start, struct point* end, struct map *m);
struct map *newMap(int height, int width) {
    struct map *m = (struct map *) malloc(sizeof(struct map));
    assert(m);
    m->height = height;
    m->width = width;

    // Note this means all values of map are 0
    int *points = (int *) calloc(height * width, sizeof(int));
    assert(points);
    m->points = (int **) malloc(width * sizeof(int *));
    assert(m->points);

    for (int i = 0; i < width; i++){
        /* Re-use sections of the memory we 
            allocated. */
        m->points[i] = points + i * height;
    }

    return m;
}

struct point *newPoint(int x, int y) {
    struct point *p = (struct point *) malloc(sizeof(struct point));
    assert(p);
    p->x = x;
    p->y = y;

    return p;
}

void freeMap(struct map *m) {
    /* We only allocate one pointer with the rest
    of the pointers in m->points pointing inside the
    allocated memory. */
    free(m->points[0]);
    free(m->points);
    free(m);
}

void printMap(struct map *m) {
    /* Print half of each row each time so we mirror the hexagonal layout. */
    int printRows = 2 * m->height;
    if(m->width < 2){
        /* If the width is less than 2, simply a vertical column so no need to print 
            last row as it will be empty. */
        printRows -= 1;
    }
    for (int i = 0; i < printRows; i++) {
        for (int j = i % 2; j < m->width; j += 2) {
            if (j == 1){
                /* For odd row, add a spacer in to place the first value after the 0th hex
                    in the printout. */
                printf("       ");
            }
            /* Default to even. Select every second column. */
            int yPos = i / 2;
            if(j % 2 != 0){
                /* If odd, numbering along height is reversed. */
                yPos = m->height - 1 - yPos;
            }
            int val = m->points[j][yPos];

            /* Print value appropriately. */
            if (val < 0){
                printf("S%.3d", val % 1000);
            } else if (val == 0){
                printf("  L  ");
            } else if (val == 100){
                printf("  A  ");
            } else {
                printf("T+%.3d ", val % 1000);
            }
            printf("          ");
        }

        printf("\n");
    }
}

/* IMPLEMENT PART A HERE */
/* Note that for the implementation in this question, you will submit
   an array that contains all the adjacent points and then add an additional
   point at the end of the array with coordinates -1, -1 to signify the
   termination of the array (similar to '\0' terminating a string) */
struct point *getAdjacentPoints(struct map *m, struct point *p) {
    struct point *ans = (struct point *) malloc(sizeof(struct point));
    assert(ans); 
    int adj_nodes = 1;
    int node_counter = 0;
    int xval = p->x;
    int yval = p->y;
     if(xval<0 || xval>=m->width || yval<0 || yval>=m->height){
        ans = realloc(ans, adj_nodes*sizeof(struct point));
        ans[node_counter].x = -1;
        ans[node_counter].y = -1;
        return ans;
    }
    for(int i = -1; i<=1;i++){
        int x_coord = p->x +i;
        for(int j = 0; j<=1;j++){
            int y_coord;
            if((i ==-1 || i==1)){
                y_coord = m->height - 1 - p->y + j;
                if(test_coord(x_coord, y_coord, m->height, m->width)){
                    adj_nodes+=1;
                    ans = realloc(ans, adj_nodes * sizeof(struct point));
                    ans[node_counter].x = x_coord;
                    ans[node_counter++].y = y_coord;
                }
            }
            else if(i == 0){
                if(j == 0){
                    y_coord = p->y + 1;
                    if(test_coord(x_coord, y_coord, m->height, m->width)){
                        adj_nodes+=1;
                        ans = realloc(ans, adj_nodes*sizeof(struct point));
                        ans[node_counter].x = x_coord;
                        ans[node_counter++].y = y_coord;
                    }
                }
                else{
                    y_coord = p->y - 1;
                    if(test_coord(x_coord, y_coord, m->height, m->width)){
                        adj_nodes+=1;
                        ans = realloc(ans, adj_nodes*sizeof(struct point));
                        ans[node_counter].x = x_coord;
                        ans[node_counter++].y = y_coord;
                }

                }
            }
        }
    }
    qsort(ans, node_counter, sizeof(struct point), comparator);
    adj_nodes+=1;
    ans = realloc(ans, adj_nodes*sizeof(struct point));
    ans[node_counter].x = -1;
    ans[node_counter].y = -1;
    return ans;
}
int comparator(const void *p, const void *q){
    int x1 = ((struct point *)p)->x;
    int x2 = ((struct point *)q)->x; 
    int y1 = ((struct point *)p)->y;
    int y2 = ((struct point *)q)->y; 
    if(x1==x2){
        return y1-y2;
    }
    else{
        return x1-x2;
    }
}
int test_coord(int x_coord, int y_coord, int height, int width){
    if(x_coord>=0 && x_coord<width && y_coord>=0 && y_coord<height){
        //printf("Test success y : %d x : %d \n", y_coord, x_coord);
        return 1;
    }
    return 0;
}
/* IMPLEMENT PART B HERE */
int mapValue(struct map *m) {
    struct point ** adj_list = getgraph(m);
    int num_points = m->height * m->width;
    int * visited = (int*)malloc(num_points* sizeof(int));
    for(int i = 0;i<num_points;i++){
        visited[i] = 0;
    }
    int total_sum = 0;
    int num_nodes[1];
    int component_sum[1];
    for(int no = 0; no<num_points;no++){
        if(visited[no] == 0 && m->points[adj_list[no][0].x][adj_list[no][0].y]>=0){
            num_nodes[0] = 0;
            component_sum[0] = 0;
            dfs(visited,adj_list, no, num_nodes, component_sum, m);
            //printf("next iteration\n");
            total_sum+= component_sum[0]*num_nodes[0]; 
        }
    }
    return total_sum;
}
void dfs(int* visited, struct point** adj_list, int no, int* num_nodes, int* component_sum, struct map *m){
    visited[no] = 1;
    //printf("(%d,%d)",adj_list[no][0].x, adj_list[no][0].y);
    //printf("(%d, %d)\n",adj_list[no][0].x,adj_list[no][0].y);
    int value = m->points[adj_list[no][0].x][adj_list[no][0].y];
    
    if(value != 100 && value>0 && component_sum[0] == 0){
        component_sum[0]+= value;
        num_nodes[0]+=1;
    }else if(value != 100 && value>0 && component_sum[0] != 0){
        component_sum[0]*=value;
        num_nodes[0]+=1;
    }
    int j = 1;
    while(adj_list [no][j].x!=-1){
        int node_num = get_node_number(adj_list[no][j].x,adj_list[no][j].y, m);
        if(visited[node_num] == 0 && get_value(adj_list[no][j].x,adj_list[no][j].y,m)>= 0){
            dfs(visited, adj_list, node_num, num_nodes, component_sum, m);
        }

        j++;
    }
}
int get_value(int x, int y, struct map *m){
    return m->points[x][y];
}
int get_node_number(int x, int y, struct map *m){
    if(x == 0){
        return y;
    }
    else{
        return (x - 0)*(m->height) + y ;
    }
}
struct point ** getgraph(struct map * m) {
    int num_points = m -> height * m -> width;
    struct point ** adj_list = malloc(sizeof(struct point *) * num_points);
    for (int i = 0; i < num_points; i++) {
        adj_list[i] = malloc(sizeof(struct point)* 8);
    }
    // Populate the first column with all possible points
    int counter = 0;
    for (int x = 0; x < m -> width ; x++) {
        for (int y = 0; y < m -> height; y++) {
            adj_list[counter][0].x = x;
            adj_list[counter][0].y = y;
            counter++;
        }
    }
    int row_counter = 0;
    for (int x_pos = 0; x_pos < m -> width; x_pos++) {
        for (int y_pos = 0; y_pos < m -> height; y_pos++) {
            struct point p;
            p.x = x_pos;
            p.y = y_pos;
            struct point * adjPoint = getAdjacentPoints(m, &p);
            //printf("reached");
            int col_counter = 1;
            //printf("(%d, %d):",p.x, p.y );
            while (adjPoint!= NULL && adjPoint ->x != -1) {
                //printf("%d\n", col_counter);
                //printf("(%d,%d)   ",adjPoint->x, adjPoint->y);
                adj_list[row_counter][col_counter].x = adjPoint -> x;
                adj_list[row_counter][col_counter].y = adjPoint -> y;
                adjPoint = adjPoint+ 1;
                col_counter+=1;
            }
            //printf("\n");
        adj_list[row_counter][col_counter].x = -1;
        adj_list[row_counter][col_counter].y = -1;
        row_counter += 1;
        }
    }
    // for(int i = 0; i< num_points; i++){
    //     int j = 0;
    //     while(adj_list[i][j].x!= -1){
    //         printf("(%d, %d)   ", adj_list[i][j].x, adj_list[i][j].y);
    //         j++;
    //     }
    //     printf("\n");
    // }
    // printf("\n\n\n");
  return adj_list;
}

/* IMPLEMENT PART D HERE */
int minTime(struct map *m, struct point *start, struct point *end) {
    struct point ** adj_list = getgraph(m);
    struct point ** mod_adj_list = add_airports(adj_list, m);
    int num_points = m->height*m->width;
    int min_time = djikstra(m, start, end, mod_adj_list);
    int time_to_airport_start = 0;
    int time_to_airport_end =0;
    struct point* closest_airport_start = get_closest_airport(start, mod_adj_list, num_points, m, &time_to_airport_start);
    struct point* closest_airport_end = get_closest_airport(end, mod_adj_list, num_points, m, &time_to_airport_end);
    int flight_time;
    if(closest_airport_start != NULL && closest_airport_end!=NULL){
        if(closest_airport_start->x ==closest_airport_end->x && closest_airport_start->y ==closest_airport_end->y ){
            flight_time = djikstra(m, closest_airport_start, end, mod_adj_list);
        }else{
            int euclid = pow(closest_airport_end->x - closest_airport_start->x, 2);
            if(15>=euclid - 85){
                flight_time = 15;
            }else{
                flight_time = euclid - 85;
            }
        }
        int total_time = time_to_airport_start + time_to_airport_end + flight_time;
        if(total_time<min_time){
            min_time = total_time;
        }
    }
    free(mod_adj_list);
    free(adj_list);
    return min_time;
}
struct point** add_airports(struct point ** adj_list, struct map *m){
    int num_points = m->height*m->width;
    struct point ** mod_adj_list = malloc(sizeof(struct point *) * num_points);
    for (int i = 0; i < num_points; i++) {
        mod_adj_list[i] = malloc(sizeof(struct point)* num_points);
    }
    //equating the two lists
    for(int i = 0; i<num_points;i++){
        int j = 0;
        while(adj_list[i][j].x!= -1){
            mod_adj_list[i][j].x = adj_list[i][j].x;
            mod_adj_list[i][j].y = adj_list[i][j].y;
            j++;
        }
        mod_adj_list[i][j].x = -1;
        mod_adj_list[i][j].y = -1;
    }
    // find all airports and add to an array
    struct point* airports = malloc(1*sizeof(struct point));
    int num_airports = 0;
    for(int i = 0; i<num_points;i++){
        if(get_value(adj_list[i][0].x,adj_list[i][0].y,m) == 100){
            airports[num_airports].x = adj_list[i][0].x;
            airports[num_airports].y = adj_list[i][0].y;
            num_airports++;
        }
    }
    for(int i = 0; i<num_points;i++){
        if(get_value(mod_adj_list[i][0].x, mod_adj_list[i][0].y, m) == 100){

            int final_index = 0;
            int append = 0;
            for(int j = 0; j<num_points;j++){
                if(mod_adj_list[i][j].x == -1){
                    final_index = j;
                    int k = 0;
                    while(k<num_airports){
                        if(!(mod_adj_list[i][0].x == airports[k].x && mod_adj_list[i][0].y == airports[k].y)){
                            mod_adj_list[i][j+append] = airports[k];
                            append++;
                        }
                        k++;
                    }
                }
            }
            mod_adj_list[i][final_index + append].x = -1;
            mod_adj_list[i][final_index + append].y = -1;
        }
    }
    free(airports);
    return mod_adj_list;

}
int djikstra(struct map *m, struct point *start, struct point *end, struct point ** adj_list){
    int num_points = m->height*m->width;
    int times[num_points];
    int min_time;
    int visited[num_points];
    for(int i = 0;i<num_points;i++){
        visited[i] = 0;
        times[i] = -1*(int)INFINITY;
    }
    times[get_node_number(start->x, start->y, m)] = 0;
    struct pq* pq = createPQ();
    insert(pq, start, 0);
    while (!isEmpty(pq)) {
        struct point* current = pull(pq);
        int curr_no = get_node_number(current->x, current->y, m);
        //printf("pulled (%d, %d)\n", current->x, current->y);
        visited[curr_no] = 1;
        if (current->x == end->x && current->y == end->y) {
            min_time = -1*times[get_node_number(end->x, end->y, m)];
            return min_time;
        }
        // Process all the neighbors of the current node
        int j = 1;
        while(adj_list [curr_no][j].x!=-1){
                int node_num = get_node_number(adj_list[curr_no][j].x,adj_list[curr_no][j].y, m);
                int time = get_time(current,&adj_list[node_num][0], adj_list, m);
                int new_time = times[curr_no] - time;
                if(!visited[node_num] && new_time > times[node_num]){
                    times[node_num] = new_time;
                    insert(pq, &adj_list[node_num][0], new_time);
                    //printf("insert (%d, %d) : %d\n", adj_list[node_num][0].x, adj_list[node_num][0].y, new_time);
                }
            j++;
        }
    }
    //printf("reached here");
    min_time = -1*times[get_node_number(end->x, end->y, m)];
    return min_time;
}
struct point* get_closest_airport(struct point *start, struct point** adj_list, int num_points, struct map *m, int* time){
    struct point* airports = malloc(num_points*sizeof(struct point));
    int num_airports = 0;
    for(int i = 0; i<num_points;i++){
        if(get_value(adj_list[i][0].x,adj_list[i][0].y,m) == 100){
            airports[num_airports].x = adj_list[i][0].x;
            airports[num_airports].y = adj_list[i][0].y;
            num_airports++;
        }
    }
    int min_airport_time = (int)INFINITY;
    int airport_number = -1;
    for(int i = 0; i<num_airports;i++){
        if(start->x != airports[i].x ||start->y != airports[i].y){
            int min_time = djikstra(m,start, &airports[i],adj_list);
            if(min_time<min_airport_time){
                min_airport_time = min_time;
                airport_number = i;
            }
        }
    }
    *time = min_airport_time;
    if(num_airports!=0){
    struct point* closest_airport = malloc(1*sizeof(struct point*));
    closest_airport->x = airports[airport_number].x;
    closest_airport->y = airports[airport_number].y;
    return closest_airport;
    }
    return NULL;

}
double get_time(struct point *start, struct point *end, struct point** adj_list, struct map *m){
    int value_start = get_value(start->x, start->y, m);
    int value_end = get_value(end->x, end->y, m);
    int time;
    if(value_start == 100 && value_end == 100){
        struct point* adjPoints = getAdjacentPoints(m, start);
        int is_adjacent = 0;
        int j = 0;
        while(adjPoints[j].x!= -1){
            if(adjPoints[j].x == start->x && adjPoints[j].y == start->y){
                is_adjacent =1;
                break;
            }
            j++;
        }
        if(is_adjacent){
            time =  5;
        }else{
            //find walking time here 
            time = (15>=pow(start->x - end->x, 2) - 85)? 15 : pow(start->x - end->x, 2) - 85;
        }
    }
    else if(value_start>0 && value_end>0){
        return 5;
    }else if(value_start < 0){
        double value_start_2 = (double) value_start;
        time = 2 + ceil((value_start_2*value_start_2)/1000);
    }
    else {
        time = 5;
    }
    return time;
}

/* IMPLEMENT PART E HERE */
int minTimeDry(struct map *m, struct point *start, struct point *end, struct point *airports, int numAirports) {
    int min_time = get_num_gaps(start, end, m);
    int time_to_airport_start = (int) INFINITY;
    int time_to_airport_end = (int)INFINITY;
    struct point closest_airport_start;
    struct point closest_airport_end;
    for(int i = 0; i <numAirports; i++){
        int time_start = get_num_gaps(start, &airports[i], m);
        int time_end = get_num_gaps(end, &airports[i], m);
        if(time_start<time_to_airport_start){
            time_to_airport_start = time_start;
            closest_airport_start.x = airports[i].x;
            closest_airport_start.y = airports[i].y;
        }
        if(time_end<time_to_airport_end){
            time_to_airport_end = time_end;
            closest_airport_end.x = airports[i].x;
            closest_airport_end.y = airports[i].y;
        }
    }
    int total_time ;
    if(numAirports>0){
        int euclid = pow(closest_airport_end.x - closest_airport_start.x, 2);
        int flight_time = (15>euclid - 85)?15:euclid - 85;
        total_time = time_to_airport_start + time_to_airport_end + flight_time;
        if(total_time<min_time){
            min_time = total_time;
        }
    }
    return min_time;
}
int get_num_gaps(struct point *start, struct point* end, struct map *m){
        int num_gaps = 0;
    while(start->x != end->x || start->y != end->y){
        int direction_x = end->x - start->x ;
        //finding absolute_direction
        int direction_y;
        if(start->x%2 ==1){
            start->y = m->height - 1 - start->y;
        }
        if(start->x%2 == 1 && end->x%2 == 1){
            direction_y = end->y - start->y;
        }
        else if(start->x%2 == 1 && end->x%2 != 1){
            int yval = m->height - 1 - start->y;
            direction_y = yval - start->y;
        }
        else if(start->x%2 != 1 && end->x%2 == 1){
            int yval = m->height - 1 - end->y;
            direction_y = yval - start->y;
        }
        else if(start->x%2 != 1 && end->x%2 != 1){
            direction_y = end->y - start->y;
        }
        if(start->x != end->x && start->y == end->y){
            num_gaps+=1;
            start->x += (direction_x>0)? 1: -1;
        }
        else if(start->x == end->x && start->y != end->y){
            num_gaps += abs(end->y) - abs(start->y);
            return 5*num_gaps;
        }
        else if(start->x != end->x && start->y != end->y){
            num_gaps+=1;
            start->x += (direction_x>0)? 1: -1;
            start->y += (direction_y>0)? 1: -1;
        }
    }
    return num_gaps*5;
}
