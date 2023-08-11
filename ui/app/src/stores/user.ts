import { defineStore } from 'pinia';
import { api } from 'src/boot/axios';
import { ref } from 'vue';

interface Endpoints {
  loginCookie: string;
  sessionLogin: string;
  sessionLogout: string;
  obtainToken?: string;
  refreshToken?: string;
  user: string;
}

interface User {
  id: number;
  email: string;
  first_name: string;
  last_name: string;
  groups: Array<string>;
}

interface ObtainTokenPayload {
  username: string;
  password: string;
}

interface UpdateTokenPayload {
  access: string;
  refresh: string;
}

export const useUserStore = defineStore('user', () => {
  const authenticated = ref<string | null>(
    localStorage.getItem('authenticated')
  );
  const user = ref<User | null>(
    JSON.parse(localStorage.getItem('user') || 'null')
  );
  const endpoints: Endpoints = {
    loginCookie: '/api/auth/cookie/',
    sessionLogin: '/api/auth/login/',
    sessionLogout: '/api/auth/logout/',
    user: '/api/auth/users/me/',
  };

  const _updateToken = (payload: UpdateTokenPayload) => {
    localStorage.setItem('authenticated', 'true');
    authenticated.value = 'true';
  };

  const _removeToken = () => {
    localStorage.removeItem('authenticated');
    authenticated.value = null;
    user.value = null;
  };

  const updateUserInfo = (payload: User) => {
    localStorage.setItem('user', JSON.stringify(payload));
    user.value = payload;
  };

  const _removeUserInfo = () => {
    localStorage.removeItem('user');
    user.value = null;
  };

  const sessionLogin = async (payload: ObtainTokenPayload) => {
    try {
      await api.get(endpoints.loginCookie);
      const resp = await api.post(endpoints.sessionLogin, payload);
      console.log('respo', resp.data);
      _updateToken(resp.data);
      await getUserInfo();
    } catch (err) {
      console.error(err);
    }
  };

  const sessionLogout = async () => {
    try {
      await removeToken();
      await api.get(endpoints.sessionLogout);
    } catch (err) {
      console.error(err);
    }
  };

  const getUserInfo = async () => {
    try {
      const resp = await api.get(endpoints.user);
      updateUserInfo(resp.data);
    } catch (err) {
      console.error(err);
    }
  };

  const removeToken = async () => {
    _removeToken();
    _removeUserInfo();
  };

  return {
    authenticated,
    user,
    sessionLogin,
    sessionLogout,
    getUserInfo,
    removeToken,
  };
});
