set backspace=indent,eol,start
set splitright
syntax on
set ignorecase
set noshowmatch
set relativenumber
set nohlsearch
set hidden
set noerrorbells
set tabstop=4 softtabstop=4
set shiftwidth=4
set expandtab
set smartindent
set nu
set wrap
set smartcase
set noswapfile
set backup
set undodir=~/.vim/undodir
set undofile
set incsearch
set termguicolors
set scrolloff=8
set nobackup

" Give more space for displaying messages.
set cmdheight=2

" Having longer updatetime (default is 4000 ms = 4 s) leads to noticeable
" delays and poor user experience.
set updatetime=50

" Don't pass messages to |ins-completion-menu|.
set shortmess+=c
set colorcolumn=80
highlight ColorColumn ctermbg=0 guibg=lightgrey
autocmd BufWritePost *.tex silent! execute "!pdflatex main.tex >/dev/null 2>&1"
set background=dark
if executable('rg')
    let g:rg_derive_root='true'
endif

au BufRead,BufNewFile *.edp setfiletype cpp
nnoremap S :%s//g<Left><Left>
nnoremap { /{<CR>
nnoremap } /}<CR>

let mapleader =" "

" Some basics:
	set nocompatible
	syntax on
	set encoding=utf-8
	"set number
	"set relativenumber

" Splits open at the bottom and right, which is non-retarded, unlike vim defaults.
	"set splitbelow
" Shortcutting split navigation, saving a keypress:
	map <C-h> <C-w>h
	map <C-j> <C-w>j
	map <C-k> <C-w>k
	map <C-l> <C-w>l
	map <leader>z ZZ
	map <leader>w :w<CR>
	map <leader>q ZQ
    nnoremap <C-f> :Files<CR>
    nnoremap <leader>pv :topleft vs<CR> :Ex<CR> :vertical resize 30<CR>
let g:netrw_browse_split = 2
    "let g:vrfr_rg = 'true'
    "let g:netrw_banner = 0
let g:netrw_winsize = 50

" Open the selected text in a split (i.e. should be a file).
	map <leader>o "oyaW:sp <C-R>o<CR>
	xnoremap <leader>o "oy<esc>:sp <C-R>o<CR>
	vnoremap <leader>o "oy<esc>:sp <C-R>o<CR>
	vnoremap <leader>t :TransSelectDirection<CR>2<CR><CR>
    map <leader>f :Goyo<CR>:set rnu<CR>

" Replace all is aliased to S.
	nnoremap S :%s//g<Left><Left>

" Open corresponding .pdf
	map <leader>p :!opout <c-r>%<CR><CR>

" Compile document
	map <leader>c :!compiler <c-r>%<CR>

"For saving view folds:
	"au BufWinLeave * mkview
	"au BufWinEnter * silent loadview

" Interpret .md files, etc. as .markdown
	let g:vimwiki_ext2syntax = {'.Rmd': 'markdown', '.rmd': 'markdown','.md': 'markdown', '.markdown': 'markdown', '.mdown': 'markdown'}

" Make calcurse notes markdown compatible:
	autocmd BufRead,BufNewFile /tmp/calcurse*,~/.calcurse/notes/* set filetype=markdown

" groff files automatically detected
	autocmd BufRead,BufNewFile *.ms,*.me,*.mom set filetype=groff

" .tex files automatically detected
	autocmd BufRead,BufNewFile *.tex set filetype=tex

" Readmes autowrap text:
	autocmd BufRead,BufNewFile *.md set tw=79

" Get line, word and character counts with F3:
	map <F3> :!wc %<CR>


" Use urlview to choose and open a url:
	:noremap <leader>u :w<Home>silent <End> !urlscan<CR>
	:noremap ,, :w<Home>silent <End> !urlscan<CR>

" Copy selected text to system clipboard (requires gvim installed):
	vnoremap <C-c> "*Y :let @+=@*<CR>
	map <C-p> "+P


	autocmd BufRead,BufNewFile /tmp/neomutt* let g:goyo_width=100

" Enable autocompletion:
	set wildmode=longest,list,full
	set wildmenu


" Automatically deletes all tralling whitespace on save.
autocmd BufWritePre * %s/\s\+$//e

" When shortcut files are updated, renew bash and ranger configs with new material:
autocmd BufWritePost ~/.key_directories,~/.key_files !bash ~/.scripts/tools/shortcuts

" Runs a script that cleans out tex build files whenever I close out of a .tex file.
autocmd VimLeave *.tex !texclear %

" Disables automatic commenting on newline:
autocmd FileType * setlocal formatoptions-=c formatoptions-=r formatoptions-=o
autocmd FileType tex autocmd FileType c,cpp,php let g:UltiSnipsExpandTrigger='<TAB>'



vmap <expr> ++ VMATH_YankAndAnalyse()
nmap ++ vip++

vnoremap K xkP`[V`]
vnoremap J xp`[V`]
vnoremap L >gv
vnoremap H <gv

set noerrorbells visualbell t_vb=
set wildmode=longest:full,full
" Use tab for trigger completion with characters ahead and navigate.
" NOTE: Use command ':verbose imap <tab>' to make sure tab is not mapped by
inoremap <silent><expr> <TAB>
      \ pumvisible() ? "\<C-n>" :
      \ <SID>check_back_space() ? "\<TAB>" :
inoremap <expr><S-TAB> pumvisible() ? "\<C-p>" : "\<C-h>"







nnoremap <Leader>b :ls<Cr>:b<Space>
nnoremap <Space>b :ls<CR>:b<Space>
inoremap <C-l> <c-g>u<Esc>[s1z=`]a<c-g>u
nmap <silent> j gj
nmap <silent> k gk
    let g:trans_directions_list = [
        \['en', 'it'],
        \['it', 'en'],
        \['', ''],
    \]
vnoremap ff :s:\\ :\\:g \| s: {:{:g \| s: _ :_:g<CR>



vnoremap <buffer> k gk
vnoremap <buffer> j gj
map <leader>gs :G<CR>
map <leader>gc :Git commit<CR>
map <leader>gp :Git push<CR>




inoremap <silent><expr> <TAB>
      \ <SID>check_back_space() ? "\<TAB>" :

function! s:check_back_space() abort
  let col = col('.') - 1
  return !col || getline('.')[col - 1]  =~# '\s'
endfunction

nnoremap <silent> K :call <SID>show_documentation()<CR>
function! s:show_documentation()
  if (index(['vim','help'], &filetype) >= 0)
    execute 'h '.expand('<cword>')
  else
    execute '!' . &keywordprg . " " . expand('<cword>')
  endif
endfunction





" Use tab for trigger completion with characters ahead and navigate.
" NOTE: Use command ':verbose imap <tab>' to make sure tab is not mapped by
inoremap <silent><expr> <TAB>
      \ pumvisible() ? "\<C-n>" :
      \ <SID>check_back_space() ? "\<TAB>" :
inoremap <expr><S-TAB> pumvisible() ? "\<C-p>" : "\<C-h>"

function! s:check_back_space() abort
  let col = col('.') - 1
  return !col || getline('.')[col - 1]  =~# '\s'
endfunction




" Source vim configuration file whenever it is saved
if has ('autocmd')          " Remain compatible with earlier versions
 augroup Reload_Vimrc       " Group name.  Always use a unique name!
    autocmd! BufWritePost .vimrc source % | echom "Reloaded " . $MYVIMRC | redraw
    autocmd! BufWritePost $MYGVIMRC if has('gui_running') | so % | echom "Reloaded " . $MYGVIMRC | endif | redraw
  augroup END
endif " has autocmd
