/**
This is written by *César Pairetti*, and is available [(here)](http://basilisk.fr/sandbox/pairetti/bag_mode/adapt_wavelet_limited.h). I keep it in the same folder for completeness. 
*/
#define TREE 1

struct Adapt_limited {
  scalar * slist1; // 第一个标量列表
  scalar * slist2; // 第二个标量列表
  double * max1;   // 第一个列表的容差
  double * max2;   // 第二个列表的容差
  int (*MLFun1)(double,double,double); // 第一个MLFun
  int (*MLFun2)(double,double,double); // 第二个MLFun
  int minlevel;    // 最小细化级别
  scalar * list;   // 需要更新的字段列表
};

trace
astats adapt_wavelet_limited (struct Adapt_limited p)
{
  scalar * listcm = NULL;

  // 合并两个slist进行restriction
  scalar * combined_slist = list_concat(p.slist1, p.slist2);

  if (is_constant(cm)) {
    if (p.list == NULL)
      p.list = all;
    restriction (combined_slist);
  }
  else {
    scalar * listr = list_concat (combined_slist, {cm});
    restriction (listr);
    free (listr);
  }
  free(combined_slist);

  astats st = {0, 0};
  scalar * listc = NULL;
  for (scalar s in p.list)
    if (!is_constant(s) && s.restriction != no_restriction)
      listc = list_add (listc, s);

  if (p.minlevel < 1)
    p.minlevel = 1;
  tree->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
  foreach_cell() {
    int cellMAX1 = p.MLFun1 ? p.MLFun1(x,y,z) : depth();
    int cellMAX2 = p.MLFun2 ? p.MLFun2(x,y,z) : depth();
    int cellMAX = max(cellMAX1, cellMAX2); // 取两者最大层级

    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
        if (cell.flags & too_coarse) {
          cell.flags &= ~too_coarse;
          refine_cell (point, listc, refined, &tree->refined);
          st.nf++;
        }
        continue;
      }
      else {
        if (cell.flags & refined) {
          cell.flags &= ~too_coarse;
          continue;
        }

        bool local = is_local(cell);
        if (!local)
          foreach_child()
            if (is_local(cell)) { local = true; break; }

        if (local) {
          static const int just_fine = 1 << (user + 3);
          
          // 处理第一个slist
          if (p.slist1) {
            int i1 = 0;
            for (scalar s in p.slist1) {
              double max1 = p.max1[i1++], sc[1 << dimension];
              int c = 0;
              foreach_child() sc[c++] = s[];
              s.prolongation (point, s);
              c = 0;
              foreach_child() {
                double e = fabs(sc[c] - s[]);
                if (e > max1 && level < cellMAX1) {
                  cell.flags &= ~too_fine;
                  cell.flags |= too_coarse;
                }
                else if ((e <= max1/1.5 || level > cellMAX1) &&
                         !(cell.flags & (too_coarse|just_fine))) {
                  if (level >= p.minlevel)
                    cell.flags |= too_fine;
                }
                else if (!(cell.flags & too_coarse)) {
                  cell.flags &= ~too_fine;
                  cell.flags |= just_fine;
                }
                s[] = sc[c++];
              }
            }
          }

          // 处理第二个slist
          if (p.slist2) {
            int i2 = 0;
            for (scalar s in p.slist2) {
              double max2 = p.max2[i2++], sc[1 << dimension];
              int c = 0;
              foreach_child() sc[c++] = s[];
              s.prolongation (point, s);
              c = 0;
              foreach_child() {
                double e = fabs(sc[c] - s[]);
                if (e > max2 && level < cellMAX2) {
                  cell.flags &= ~too_fine;
                  cell.flags |= too_coarse;
                }
                else if ((e <= max2/1.5 || level > cellMAX2) &&
                         !(cell.flags & (too_coarse|just_fine))) {
                  if (level >= p.minlevel)
                    cell.flags |= too_fine;
                }
                else if (!(cell.flags & too_coarse)) {
                  cell.flags &= ~too_fine;
                  cell.flags |= just_fine;
                }
                s[] = sc[c++];
              }
            }
          }

          foreach_child() {
            cell.flags &= ~just_fine;
            if (!is_leaf(cell)) {
              cell.flags &= ~too_coarse;
              if (level >= cellMAX)
                cell.flags |= too_fine;
            }
            else if (!is_active(cell))
              cell.flags &= ~too_coarse;
          }
        }
      }
    }
    else continue;
  }
  mpi_boundary_refine (listc);

  for (int l = depth(); l >= p.minlevel; l--) {
    foreach_cell()
      if (!is_boundary(cell)) {
        if (level == l) {
          if (!is_leaf(cell)) {
            if (cell.flags & refined)
              cell.flags &= ~(refined|too_fine);
            else if (cell.flags & too_fine) {
              if (is_local(cell) && coarsen_cell (point, listc))
                st.nc++;
              cell.flags &= ~too_fine;
            }
          }
          if (cell.flags & too_fine)
            cell.flags &= ~too_fine;
          else if (aparent(0).flags & too_fine)
            aparent(0).flags &= ~too_fine;
          continue;
        }
        else if (is_leaf(cell))
          continue;
      }
    mpi_boundary_coarsen (l, too_fine);
  }
  free (listc);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (p.list);
  free (listcm);

  return st;
}