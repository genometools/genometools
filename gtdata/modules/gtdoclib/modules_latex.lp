\section*{Table of Modules}
\begin{itemize}
<%for _, module in ipairs(modules) do%>
  \item Module \texttt{<%=module%>}\dotfill page \pageref{<%=module%>}
<%end%>
\end{itemize}
