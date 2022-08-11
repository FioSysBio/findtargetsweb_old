from django.conf.urls import url
from . import views

from django.conf import settings
from django.conf.urls.static import static

app_name = "FindTargetsWEB"

urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'passo1/$', views.passo1, name='passo1'),
    url(r'passo2/$', views.passo2, name='passo2'),
    url(r'download/$', views.download, name='download'),
]+ static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
